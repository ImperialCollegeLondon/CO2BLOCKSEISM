% read inout parameters
% Note that seismicity may occur in formations overlying or underlying and 
% in hydraulic conncetion with the target storage layer. Thus, depth,
% pressure, stress and fault properties are entered separately for the
% seismic layer. 
function [thick,area_res,perm,por,dens_c,visc_c,visc_w,compr,rc,...
    gamma,delta,omega,depth_seismic,pres_grad_seismic,Sv_grad_seismic,...
    Shmin_grad_seismic,SHmax_grad_seismic,SHmax_dir_seismic,fault_friction]...
    = read_data(path,name)
 
    % -- default parameters (used in case they are not provided)
    def_litho_grad = 23 ;               % lithostatic gradient [MPa/km]
    def_hydro_grad = 10 ;               % hydrostatic gradient [MPa/km]
    def_temp_grad = 33 ;                % temperaturegradient [C/km]
    def_temp_surface = 15 ;             % surface temperature [C]
    def_k0 = 0.6 ;                      % default stress ratio Sh/Sv  
    def_rock_friction_angle = 30 ;      % default rock friction angle [deg]
    def_rock_cohesion = 0;              % default rock cohesion [MPa]
    def_cr = 5*10^-4 ;                  % default rock compressibility [MPa^-1]
    def_cw = 3*10^-4 ;                  % default water compressibility [MPa^-1]
    def_salinity = 180000 ;             % default salinity [ppm]
    def_fault_friction = 0.58;          % default fault friction coefficient
    def_SHmax_dir_seismic = 90;         % default direction for SHmax in the seismic layer [deg]
    % --

    % -- read data
    data = readtable([path,'\',name],'Sheet','Model parameters','Range','A1:AD2'); 

    % General data used for calculating injectivity in "calculate" 
    domain_type = char(data{1,2});                      % domain confinement; 'open' or 'closed'
    depth = double(data{1,3});                          % shallowest depth of reservoir [m]
    depth_mean = double(data{1,4});                     % mean depth of reservoir [m]
    thick = double(data{1,5});                          % thickness of reservoir [m]
    area_res = double(data{1,6});                       % area of reservoir [km^2]
    perm = double(data{1,7})*10^-15;                    % intrinsic permeability [m^2]
    por = double(data{1,8});                            % porosity [-]
    cr = double(data{1,9})/1e6;                         % rock compressibility [1/Pa]
    cw = double(data{1,10})/1e6;                        % water compressibility [1/Pa]
    dens_c = double(data{1,11})*1e3;                    % Density of CO2 [kg/m^3]
    visc_c = double(data{1,12})/1e3;                    % Viscosity of CO2 [Pa.s] 
    visc_w = double(data{1,13})/1e3;                    % Viscosity of water[Pa.s] 
    pres_grad = double(data{1,14})/1e3;                 % pressure gradient [MPa/m]
    temp_surf = double(data{1,15});                     % surface temperature [C]
    temp_grad = double(data{1,16})/1e3;                 % temperature gradient [C/m]
    Sv_grad = double(data{1,17})/1e3;                   % overburden stress gradient (MPa/m)
    Shmin_grad = double(data{1,18})/1e3;                % Shmin stress gradient (MPa/m)
    SHmax_grad = double(data{1,19})/1e3;                % SHmax stress gradient (MPa/m)
    salinity = double(data{1,20})/1e6;                  % aquifer salinity [ppm/1e6]
    rock_friction = double(data{1,21});                 % rock friction angle  [deg]
    rock_cohesion = double(data{1,22});                 % rock cohesion coefficient [MPa]
    rock_tens_strength = double(data{1,23});            % rock tensile strength [MPa]

    % Seismicity data
    depth_seismic = double(data{1,24});                 % depth at which seismicity may occur [m]
    pres_grad_seismic = double(data{1,25})/1e3;         % pore pressure gradient (MPa/m)
    Sv_grad_seismic = double(data{1,26})/1e3;           % overburden stress gradient (MPa/m)
    Shmin_grad_seismic = double(data{1,27})/1e3;        % Shmin stress gradient (MPa/m)
    SHmax_grad_seismic = double(data{1,28})/1e3;        % SHmax stress gradient (MPa/m)
    SHmax_dir_seismic = double(data{1,29});             % SHmax stress direction (deg)
    fault_friction = double(data{1,30});                % fault friction coefficient
    % --

    switch domain_type
        % rc is the equivalent circle-shape reservoir radius
        case 'Open'
            domain_type = 'open';
            rc = inf;
        case 'open'
            domain_type = 'open';
            rc = inf;
        case 'Closed'
            domain_type = 'closed';
            rc = sqrt(area_res*10^6/pi);    
        case 'closed'
            domain_type = 'closed';
            rc = sqrt(area_res*10^6/pi);  
    end  

    %%% calculate some parameters if not given
    if depth == 0 ||  isnan(depth)
        depth = depth_mean - thick/2 ;
    end

    if pres_grad == 0 ||  isnan(pres_grad)
        pres_grad = def_hydro_grad/1e3 ;
    end
    
    if temp_surf == 0 ||  isnan(temp_surf)
        temp_surf = def_temp_surf ;
    end

    if temp_grad == 0 ||  isnan(temp_grad)
        temp_grad = def_temp_grad/1e3 ;
    end

    if Sv_grad == 0  ||  isnan(Sv_grad)
        Sv_grad = def_litho_grad/1e3;
    end

    if Shmin_grad == 0  ||  isnan(Shmin_grad) && SHmax_grad>0
        Shmin_grad = SHmax_grad;
    end

    if SHmax_grad == 0  ||  isnan(SHmax_grad) && Shmin_grad>0
        SHmax_grad = Shmin_grad;
    end

    if SHmax_grad == 0  ||  isnan(SHmax_grad) && Shmin_grad == 0  ||  isnan(Shmin_grad)
        Shmin_grad = def_k0*Sv_grad ;
        SHmax_grad = def_k0*Sv_grad ;
    end

    if rock_friction == 0 ||  isnan(rock_friction)
        rock_friction = def_rock_friction_angle ;
    end
    
    if rock_cohesion == 0 ||  isnan(rock_cohesion)
        rock_cohesion = def_rock_cohesion ;
    end

    if rock_tens_strength== 0 ||  isnan(rock_tens_strength)
        rock_tens_strength = rock_cohesion/2;
    end
                
    if cr == 0 ||  isnan(cr)
        cr = def_cr/1e6;
    end

    if cw == 0 ||  isnan(cw)
        cw = def_cw/1e6;
    end

    if salinity == 0 || isnan(salinity)
        salinity = def_salinity/1e6 ;
    end

    % Seismicity assessment data
    if depth_seismic == 0 || isnan(depth_seismic)
        depth_seismic = depth_mean + thick/2 ;
    end

    if pres_grad_seismic == 0 || isnan(pres_grad_seismic)
        pres_grad_seismic = pres_grad ;
    end

    if Sv_grad_seismic == 0 || isnan(Sv_grad_seismic)
        Sv_grad_seismic = Sv_grad ;
    end
    
    if Shmin_grad_seismic == 0 || isnan(Shmin_grad_seismic)
        Shmin_grad_seismic = Shmin_grad ;
    end

    if SHmax_grad_seismic == 0 || isnan(SHmax_grad_seismic)
        SHmax_grad_seismic = SHmax_grad ;
    end
    
    if SHmax_dir_seismic == 0 || isnan(SHmax_dir_seismic)
        SHmax_dir_seismic = def_SHmax_dir_seismic ;
    end

    if fault_friction == 0 || isnan(fault_friction)
        fault_friction = def_fault_friction ;
    end

    %---------------------------------------------------------------------
        
    depth_bottom = depth_mean + thick/2 ;                % bottom reservoir depth [m]
    pres_mean = pres_grad*depth_mean ;                   % pressure at mid reservoir depth [MPa]
    T_mean = temp_surf + temp_grad*depth_mean;          % temperature at mid reservoir depth [C]

    % density and viscosity calculation at mean reservoir p&T
    if dens_c == 0 ||  isnan(dens_c)
        [~,dens_c,~] = eos(T_mean, pres_mean,salinity,0);
    end

    if visc_c == 0 || isnan(visc_c)
        [~,~,visc_c] = eos(T_mean, pres_mean,salinity, dens_c);
    end

    if visc_w == 0 || isnan(visc_w) 
        [visc_w,~,~] = eos(T_mean, pres_mean, salinity,0) ; 
    end

    %%%%%%%%% calculate some useful parameters 
    gamma = visc_c/(visc_w);                                                % Non-dimensional viscosity ratio[-]
    delta = (visc_w-visc_c)/visc_w;  
    omega = (visc_c+visc_w)/(visc_c-visc_w)*log(sqrt(visc_c/visc_w))-1 ;
    compr = cr+por*cw ;                                                     %total compressibility  [1/Pa]
end
