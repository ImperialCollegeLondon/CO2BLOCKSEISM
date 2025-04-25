function [Mesh_grid,d_list,d_max_square,well_list,p_fault,p_2Dgrid] = ...
    hydrogeology(fpath,fname,model_width,model_height,...
    fault_coord_x,fault_coord_y,dist_min,dist_max,nr_dist,rw,nr_well_max,dt,...
    time_project,Q_M_each,nr_fault,n_2Dplot_grids_x,n_2Dplot_grids_y)

    %read data
    [thick,area_res,perm,por,dens_c,visc_c,visc_w,compr,rc,gamma,~,omega,~,~,~,~,~,~,~,~,~,~,~,~] = read_data(fpath,fname); %#ok<*LLMNC>

    % initialize
    if strcmp(nr_well_max,'auto')                                           % calculate maximum well number if not set
        nr_well_max = floor(area_res/(dist_min^2));
    end

    if strcmp(dist_max,'auto')                                              % calculate maximum interwell distance if not set
        dist_max = sqrt(2*area_res)/2;                                      % radius of the circle circumscribing square-shape reservoir area
    end

    d_list = linspace(dist_min,dist_max,nr_dist) ;                          % inter-well distance list

    well_list = [];
    d_max_square = [];
    w_id = 0 ;

    dimension_x = 1000*model_width;                      % Model dimension in x direction [m]
    dimension_y = 1000*model_height;                     % Model dimension in y direction [m]

    % calculate grid and distances
    for x_grid_num = 1:sqrt(nr_well_max)                                    % number of wells on a horizontal row
        if x_grid_num*(x_grid_num+1) < nr_well_max
            plus = 1;
        else
            plus = 0;
        end
        for y_grid_num = x_grid_num:x_grid_num+plus                         % number of wells on a  vertical row
            w_id = w_id +1 ;                                                % well number scenario (configuration) ID
            w = x_grid_num*y_grid_num ;                                     % well number for each scenario
            well_list(1,w_id) = w;                                          % store in vector
            well_list(2,w_id) = x_grid_num;                                 % number of x grids stored in vector
            well_list(3,w_id) = y_grid_num;                                 % number of y grids stored in vector
            d_max_square(w_id) = sqrt(area_res)/y_grid_num;                 % maximum interwell distance for each well number [km]

            % grid coordinates for contour plots of overpressure distribution
            x_grid = linspace(0,dimension_x,n_2Dplot_grids_x);
            y_grid = linspace(0,dimension_y,n_2Dplot_grids_y);
            [X,Y] = meshgrid(x_grid,y_grid);
            Mesh_grid{w_id}(:,:,1) = X;
            Mesh_grid{w_id}(:,:,2) = Y;

            for d = 1:nr_dist
                if d_list(d) <= d_max_square(w_id)        % calculate pressure distribution only for allowable interwell distances
                    distance= d_list(d)*1000;

                    % (x,y) coordinate of each wellbore [m]

                    wells_coord_x = (model_width - sqrt(area_res))*1000/2 + repmat((0:distance:distance*x_grid_num-1) + (y_grid_num*1000*d_max_square(w_id)-(x_grid_num-1)*distance)/2,y_grid_num,1);
                    wells_coord_y = (model_height - sqrt(area_res))*1000/2 + repmat(transpose((distance*(y_grid_num-1):-distance:0) + (y_grid_num*1000*d_max_square(w_id)-(y_grid_num-1)*distance)/2),1,x_grid_num);

    %%%%%%%%%%%%%%%%% Calculate pressure build-up at each fault using the
                    % maximum allowable injection rate for each scenario (configuration, distance)
                    % calculated from the calculate module
                    % Q_v = Q_M_each*1e9/dens_c/365/86400;                   % injection rate per well [m3/s]

                    % Calculations of pressure buildup on faults
                    for fault_number = 1:nr_fault
                        dist_vec_x = wells_coord_x - fault_coord_x(fault_number);              % distance in x between the selected fault and wellbores [m]
                        dist_vec_y = wells_coord_y - fault_coord_y(fault_number);              % distance in y between the selected fault and wellbores [m]
                        dist_vec   = sqrt(dist_vec_x.^2+dist_vec_y.^2) ;                             % distance between the selected fault and all wells [m]

                        time_count = 0;
                        for time_step = dt:dt:time_project
                            time_count = time_count + 1;
                            time = time_step*86400*365 ;                               % injection time [sec]

                            p_sup = 0;
                            for i = 1:w  % w is the number of the wells in each scenario
                                r = dist_vec(i);
                                if r<rw
                                    r=rw;
                                end

                                % finding where t in the Q_M timing list is
                                % located
                                bigger_value_location = find((Q_M_each(1,:)-time_step)>=0);
                                if bigger_value_location>0
                                    previous_rate_no = bigger_value_location(1)-1;
                                else
                                    [~,previous_rate_no] = size(Q_M_each);
                                end

                                % calculating CO2 volume injected so far &
                                % csi for Q0
                                vol_CO2 = time_step*Q_M_each(2,1)*1e9/dens_c;

                                csi = sqrt(vol_CO2/pi/por/thick);                          % average plume extension [m]
                                psi = exp(omega)*csi;                                      % equivalent plume extension [m]

                                R_influence = sqrt(2.246*perm*time/(visc_w*compr));        % pressure propagation radius for the time of injection [m]
                                p_c = (Q_M_each(2,1)*1e9/dens_c/365/86400*visc_w)/(2*pi*thick*perm)/1e6;        % characteristic pressure [MPa]

                                % Baseline pressure increase by Q0
                                Delta_p = Analytical_solution(r,R_influence,psi,rc,gamma)*p_c;    % overpressure according to Nordbotten and Celia solution [MPa]

                                if previous_rate_no > 1
                                    for v= 2:previous_rate_no
                                        vol_CO2 = vol_CO2 + (time_step-Q_M_each(1,v))*(Q_M_each(2,v)-Q_M_each(2,v-1))*1e9/dens_c;  % volume of CO2 injected
                                        csi = sqrt(vol_CO2/pi/por/thick);                          % average plume extension [m]
                                        psi = exp(omega)*csi;                                      % equivalent plume extension [m]
                                        p_c = ((Q_M_each(2,v) - Q_M_each(2,v-1))*1e9/dens_c/365/86400*visc_w)/(2*pi*thick*perm)/1e6;         % characteristic pressure [MPa]
                                        R_influence = sqrt(2.246*perm*(time_step - Q_M_each(1,v))*86400*365/(visc_w*compr));       % pressure propagation radius for the time of injection [m]

                                        % the case where the CO2 plume from previous injection rates has propagated to a large distance
                                        % in this case, R is recalculated by mu_co2 and pressure is calculated by the function gamma*log(R_new/r)
                                        if R_influence <= psi
                                            R_influence = sqrt(2.246*perm*(time_step - Q_M_each(1,v))*86400*365/(visc_c*compr));
                                            Delta_p = Delta_p + Analytical_solution(r,R_influence,R_influence,rc,gamma)*p_c;    % overpressure according to Nordbotten and Celia solution [MPa]
                                        else
                                            Delta_p = Delta_p + Analytical_solution(r,R_influence,psi,rc,gamma)*p_c;    % overpressure according to Nordbotten and Celia solution [MPa]
                                        end
                                    end
                                end
                                p_sup = p_sup +  Delta_p ;                             % superposed overpressure for the center wellbore [MPa]
                            end

                            p_fault{w_id}{d}(fault_number,time_count) = p_sup ;        % pressure buildup on each fault for particular distance and well configuration [MPa]
                        end
                    end


    %%%%%%%%%%%%%%% Overpressure calculations over the reservoir domain
                    for x = 1:n_2Dplot_grids_y         % x is the row number
                        for y = 1:n_2Dplot_grids_x     % y is the column number
                            dist_vec_x = wells_coord_x - X(x,y);                       % distance in x between the selected grid point and wellbores [m]
                            dist_vec_y = wells_coord_y - Y(x,y);                       % distance in y between the selected grid point and wellbores [m]
                            dist_vec   = sqrt(dist_vec_x.^2 + dist_vec_y.^2) ;         % distance between the selected grid point and all wells [m]

                            time_count = 0;
                            for time_step = dt:dt:time_project
                                time_count = time_count + 1;
                                time = time_step*86400*365 ;                               % injection time [sec]

                                p_sup = 0;
                                for i = 1:w
                                    r = dist_vec(i);
                                    if r<rw
                                        r=rw;
                                    end

                                    % finding where t in the Q_M timing list is
                                    % located
                                    bigger_value_location = find((Q_M_each(1,:)-time_step)>=0);
                                    if bigger_value_location>0
                                        previous_rate_no = bigger_value_location(1)-1;
                                    else
                                        [~,previous_rate_no] = size(Q_M_each);
                                    end

                                    % calculating CO2 volume injected so far & csi
                                    vol_CO2 = time_step*Q_M_each(2,1)*1e9/dens_c;              % m3

                                    csi = sqrt(vol_CO2/pi/por/thick);                          % average plume extension [m]
                                    psi = exp(omega)*csi;                                      % equivalent plume extension [m]
                                    R_influence = sqrt(2.246*perm*time/(visc_w*compr));        % pressure propagation radius for the time of injection [m]
                                    p_c = (Q_M_each(2,1)*1e9/dens_c/365/86400*visc_w)/(2*pi*thick*perm)/1e6;        % characteristic pressure [MPa]

                                    % Baseline pressure increase by Q0
                                    Delta_p = Analytical_solution(r,R_influence,psi,rc,gamma)*p_c;       % overpressure according to Nordbotten and Celia solution [MPa]

                                    if previous_rate_no > 1
                                        for v= 2:previous_rate_no
                                            vol_CO2 = vol_CO2 + (time_step-Q_M_each(1,v))*(Q_M_each(2,v)-Q_M_each(2,v-1))*1e9/dens_c;   % CO2 volume for injection until the time step
                                            csi = sqrt(vol_CO2/pi/por/thick);                          % average plume extension [m]
                                            psi = exp(omega)*csi;                                      % equivalent plume extension [m]
                                            p_c = ((Q_M_each(2,v) - Q_M_each(2,v-1))*1e9/dens_c/365/86400*visc_w)/(2*pi*thick*perm)/1e6;         % characteristic pressure [MPa]
                                            R_influence = sqrt(2.246*perm*(time_step - Q_M_each(1,v))*86400*365/(visc_w*compr));       % pressure propagation radius for the time of injection [m]

                                            % the case where the CO2 plume from previous injection rates has propagated to a large distance
                                            % in this case, R is recalculated by mu_co2 and pressure is calculated by the function gamma*log(R_new/r)
                                            if R_influence <= psi
                                                R_influence = sqrt(2.246*perm*(time_step - Q_M_each(1,v))*86400*365/(visc_c*compr));
                                                Delta_p = Delta_p + Analytical_solution(r,R_influence,R_influence,rc,gamma)*p_c;
                                            else
                                                Delta_p = Delta_p + Analytical_solution(r,R_influence,psi,rc,gamma)*p_c;    % overpressure according to Nordbotten and Celia solution [MPa]
                                            end
                                        end
                                    end
                                    p_sup = p_sup +  Delta_p ;                             % superposed overpressure for the center wellbore [MPa]
                                end

                                p_2Dgrid{w_id}{d}{time_count}(x,y) = p_sup ;               % pressure buildup across the reservoir for different times, a particular well configuration and distance [MPa]
                            end
                        end
                    end
                end
            end
        end
    end
end
