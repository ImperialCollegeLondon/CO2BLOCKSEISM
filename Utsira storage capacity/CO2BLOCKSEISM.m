% This is a preliminary version of the tool "CO2BLOCKSEISM", which provides
% estimate of the CO2 storage capacity of a geological reservoir constrained
% by the potential for induced seismicity. Injection sites are distributed
% on rectangular grids at distance d in either m×m or ‎m×(m+1) configurations. ‎
% The tool is here applied to the screening of storage capacity of the
% Utsira Formation in the North Sea

% Please, read the file README for more information about the tool.
% This software is free. Please cite CO2BLOCKSEISM as:
% https://github.com/imanrahimzadeh/CO2BLOCKSEISM
%
% Kivi, I.R., De Simone, S. and Krevor, S., 2024. A simplified physics model
% for estimating subsurface CO2 storage resources constrained ‎by % fault slip
% potential. Preprint, EarthArXiv. DOI: https://doi.org/10.31223/X5KM65
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIXME: deduplicate this script by introducing functions

clearvars; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- data file directory and name
fpath = pwd;                        % directory of the input data file
addpath(fpath)

fname = 'Input-Utsira.xlsx';               % name of the input data file

Fault_data = xlsread("Faults.xlsx");                            % Fault attributes

ref_X = 46850;                       % X coordinate of the reference point [m]
ref_Y = 6443400;                     % Y coordinate of the reference point [m]

model_width = 120;                  % model dimension in the longitude direction [km]
model_height = 120;                 % model dimension in the latitude direction [km]

% General wellbore grid parameters
dist_min = 5 ;                      % minimum inter-well distance [km]; default= 5 km
dist_max = 'auto';                  % maximum inter-well distance [km]. Set a number or 'auto' if you prefer automatic calculation
nr_dist = 15;                       % number of inter-well distances to explore; default= 15
nr_well_max ='auto';                % maximum number of wells. Set a number or 'auto' if you prefer automatic calculation
rw = 0.1 ;                          % well radius [m]

time_project = 50 ;                 % projects duration [years]

% Assigning parameters for the probabilistic analysis of seismicity
% assigning distribution functions ('Nor': normal, 'Uni': Uniform)
f_dist_p = 'Uni';
f_dist_Sv = 'Nor';
f_dist_Shmin = 'Nor';
f_dist_SHmax = 'Nor';
f_dist_dir = 'Nor';
f_dist_dip = 'Nor';
f_dist_azi = 'Nor';
f_dist_mu = 'Nor';

% Distribution of uncertain geomechanical parameters
% Option 1: Uniform distribution function
% +/-a around the average value
% or
% Option 2: Normal distribution function
% generating values within +-3a of the average value (a: standard deviation)
a_p_grad = 0.2 ;                    % pressure gradient [MPa/km]
a_Sv_grad = 0.2 ;                   % vertical stress gradient [MPa/km]
a_Shmin_grad = 1 ;                  % Shmin stress gradient [MPa/km]
a_SHmax_grad = 0.1 ;                  % SHmax stress gradient [MPa/km]   0.5
a_SHmax_dir = 10 ;                  % Azimuth of SHmax [deg]
a_fault_dip = 18 ;                   % Fault dip [deg]
a_fault_azi = 2 ;                   % Fault azimuth [deg]
a_fault_mu = 0.03;                  % Fault friction coefficien

n_MC = 5000 ;                        % Number of Monte-Carlo Simulations

% Assigning fault attributes
% Dip and Azimuth
fault_dip = Fault_data(:,5);
fault_azi = Fault_data(:,4);

% fault length [m]
fault_length = Fault_data(:,3);

% fault coordinate (centroid) [m]
fault_coord_x = Fault_data(:,1) - ref_X;
fault_coord_y = Fault_data(:,2) - ref_Y;

nr_fault = length(fault_dip);


% Assigning spatial and temporal domain resolutions
dt= 1;                            % time steps overwhich pressure calculations are made [y]; t: 1:dt:time_project
n_2Dplot_grids_x = 201;             % number of grid points along x for pressure calculation across the reservoir
n_2Dplot_grids_y = 201;             % number of grid points along y for pressure calculation across the reservoir

% Assigning the injection scheme

Q_M = [0 ;...         % Injection rate per well [Mt/y]
          10];        % Different number of injection plans can be entered
                      % Each injection plan has two rows: time and rates, respectively
                      % the time corresponds to the starts of using each rate

%%%%%%%%%%%%%%%%%%%%%%%%%% END OF INPUT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%    CALCULATIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[deltap_cr,deltap_cr_ref,r_p,r_Sv,r_Shmin,r_SHmax,r_SHmax_dir,r_fault_mu,r_fault_dip,r_fault_azi,...
    pd_p,pd_Sv,pd_Shmin,pd_SHmax,pd_SHmax_dir,pd_fault_mu,pd_fault_dip,pd_fault_azi]...
    = Geomech_prob(fpath,fname,fault_dip,fault_azi,f_dist_p,f_dist_Sv,f_dist_Shmin,f_dist_SHmax,...
    f_dist_dir,f_dist_dip,f_dist_azi,f_dist_mu,a_p_grad,a_Sv_grad,a_Shmin_grad,a_SHmax_grad,...
    a_SHmax_dir,a_fault_dip,a_fault_azi,a_fault_mu,nr_fault,n_MC);
disp('probabilistic assessment done')


[Mesh_grid,d_list,d_max_square,well_list,p_fault,p_2Dgrid] = hydrogeology(fpath,fname,...
    model_width,model_height,fault_coord_x,fault_coord_y,dist_min,dist_max,nr_dist,...
    rw,nr_well_max,dt,time_project,Q_M,nr_fault,n_2Dplot_grids_x,n_2Dplot_grids_y);
disp('hydrogeology calculations done')

%%%%%%%%%%%%%%%%%%%%%%%%% END OF CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
