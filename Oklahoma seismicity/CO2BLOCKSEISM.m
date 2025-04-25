% This is a preliminary version of the tool "CO2BLOCKSEISM", which provides
% estimate of the CO2 storage capacity of a geological reservoir constrained
% by the potential for induced seismicity. The tool is demonstrated here
% by modeling induced seismicity as a result of wastewater disposal in
% Oklahoma, US.

% Please, read the file README for more information about the tool.
% This software is free. Please cite CO2BLOCKSEISM as:
% https://github.com/imanrahimzadeh/CO2BLOCKSEISM
%
% Kivi, I.R., De Simone, S. and Krevor, S., 2024. A simplified physics model
% for estimating subsurface CO2 storage resources constrained ‎by % fault slip
% potential. Preprint, EarthArXiv. DOI: https://doi.org/10.31223/X5KM65
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%#ok<*LLMNC>

clearvars; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- data file directory and name
fpath = pwd;                        % directory of the input data file
addpath(fpath)

fname = 'Input_OK.xlsx';               % name of the input data file

well_data = readmatrix("Inj_rate.xlsx");              % Monthly injection rate data of the wellbores
Fault_data = readmatrix("Faults.xlsx");                  % Fault attributes

Study_area = [35.2 -99.4; 35.2 -96.15;...
    37.5 -96.15; 37.5 -99.4;35.2 -99.4];          % Coordinate of the study area


ref_lat = 35;                       % latitude of the reference (bottom-left corner) [º]
ref_lon = -99.5;                    % longitude of the reference (bottom_left corner) [º]

model_width = 300;                  % model dimension in the longitude direction [km]
model_height = 300;                 % model dimension in the latitude direction [km]

%-- setting parameters
% Injectivity parameters
rw = 0.1 ;                          % well radius [m]

time_project = 18 ;                 % projects duration [years]

%%%%%%%%%%%%%%%%   Probabilistic geomechanics parameters  %%%%%%%%%%%%%%%%%
% Assigning parameters for the probabilistic analysis of fault slip
stress_criticality = 'Y';         % 'Y': critically stressed
                                  % (only for SS regime: SHmax is calculated from other parameters)
                                  % 'N': all stress components need to be specified

dp_crit = 0;                      % pressure required to reactivate an optimally oriented fault [MPa]
                                  % only considered for stress_criticality = 'Y'

% assigning distribution functions ('Nor': normal, 'Uni': Uniform)
f_dist_p = 'Uni';
f_dist_Sv = 'Nor';
f_dist_Shmin = 'Nor';
f_dist_SHmax = 'Nor';
f_dist_dir = 'Nor';
f_dist_dip = 'Nor';
f_dist_azi = 'Nor';
f_dist_mu = 'Nor';

% Uniform distribution function
% +/-a around the average value (already read from the Input data file or below)
% or
% Normal distribution function
% generating values within +-3a of the average value (a: standard deviation)
a_p_grad = 0.6 ;                    % pressure gradient [MPa/km]
a_Sv_grad = 0.5 ;                   % vertical stress gradient [MPa/km]
a_Shmin_grad = 1 ;                  % Shmin stress gradient [MPa/km]
a_SHmax_grad = 1.5 ;                % SHmax stress gradient [MPa/km]   2.152 by walsh
a_SHmax_dir = 1.5 ;                   % Azimuth of SHmax [deg]
a_fault_dip = 5 ;                   % Fault dip [deg]
a_fault_azi = 2 ;                   % Fault azimuth [deg]
a_fault_mu = 0.026;                 % Fault friction coefficien

n_MC = 5000 ;                        % Number of Monte-Carlo Simulations

% Assigning fault attributes (they are known in this case)
fault_dip = Fault_data(:,5);
fault_azi = Fault_data(:,4);

% fault length [m]
fault_length = Fault_data(:,3);

% fault coordinate (centroid) [m]
fault_coord_x = Fault_data(:,1);
fault_coord_y = Fault_data(:,2);

nr_fault = length(fault_dip);

% Assigning spatial and temporal domain resolutions
dt= 1/12;                          % time steps overwhich pressure calculations are made [y];
                                   % fixed at one month for this problem as the input parameters are all in monthly intervals
n_2Dplot_grids_x = 101;             % number of grid points along x for pressure calculation across the reservoir
n_2Dplot_grids_y = 101;             % number of grid points along y for pressure calculation across the reservoir

% Loading injection rate data [bbl/month]
Inj_data = well_data(:,5:220);

% Cumulative injected volume [m3]
cum_water = well_data(:,4)*0.159;

% well coordinates [m]
wells_coord_x = 111.320*1000*(well_data(:,3) - ref_lon).*cos(deg2rad(well_data(:,2)));
wells_coord_y = 110.574*1000*(well_data(:,2) - ref_lat);

%%%%%%%%%%%%%%%%%%%%%%%%%% END OF INPUT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
% calculating slip probability curves
[deltap_cr,deltap_cr_ref,r_p,r_Sv,r_Shmin,r_SHmax,r_SHmax_dir,r_fault_mu,r_fault_dip,r_fault_azi,...
    pd_p,pd_Sv,pd_Shmin,pd_SHmax,pd_SHmax_dir,pd_fault_mu,pd_fault_dip,pd_fault_azi]...
    = Geomech_prob(fpath,fname,stress_criticality,dp_crit,fault_dip,fault_azi,f_dist_p,...
    f_dist_Sv,f_dist_Shmin,f_dist_SHmax,f_dist_dir,f_dist_dip,f_dist_azi,f_dist_mu,...
    a_p_grad,a_Sv_grad,a_Shmin_grad,a_SHmax_grad,a_SHmax_dir,a_fault_dip,a_fault_azi,...
    a_fault_mu,nr_fault,n_MC);
disp('probabilistic assessment done')

% calculating injection overpressure
[Mesh_grid,p_fault,p_2Dgrid] = hydrogeology(fpath,fname,model_width,model_height,wells_coord_x,wells_coord_y,...
    fault_coord_x,fault_coord_y,rw,Inj_data,dt,time_project,nr_fault,n_2Dplot_grids_x,n_2Dplot_grids_y);
disp('hydrogeology calculations done')

%%%%%%%%%%%%%%%%%%%%%%%%% END OF CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
