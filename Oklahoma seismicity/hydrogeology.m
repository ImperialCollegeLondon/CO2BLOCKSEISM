%#ok<*LLMNC>
function [Mesh_grid,p_fault_reservoir,p_2Dgrid] = ...
    hydrogeology(fpath,fname,model_width,model_height,wells_coord_x,wells_coord_y,...
    fault_coord_x,fault_coord_y,rw,Q_M,dt,time_project,nr_fault,...
    n_2Dplot_grids_x,n_2Dplot_grids_y)

    % read data
    [thick,~,perm,~,~,~,visc_w,compr,rc,~,~, ~,~,~,~,~,~,~,~] = read_data(fpath,fname);


    % grid coordinates for contour plots of overpressure distribution
    x_grid = linspace(0,(model_width*1000),n_2Dplot_grids_x);
    y_grid = linspace(0,(model_height*1000),n_2Dplot_grids_y);
    [X,Y] = meshgrid(x_grid,y_grid);
    Mesh_grid(:,:,1) = X;
    Mesh_grid(:,:,2) = Y;

    % Calculations of pressure buildup on faults
    for fault_number = 1:nr_fault
        dist_vec_x = wells_coord_x - fault_coord_x(fault_number);              % distance in x between the selected fault and wellbores [m]
        dist_vec_y = wells_coord_y - fault_coord_y(fault_number);              % distance in y between the selected fault and wellbores [m]
        dist_vec   = sqrt(dist_vec_x.^2+dist_vec_y.^2) ;                       % distance between the selected fault and all wells [m]

        time_count = 0;
        for time_step = dt:dt:time_project
            time_count = time_count + 1;
            time = time_step*86400*365 ;                               % injection time [sec]

            p_sup = 0;
            for i = 1:length(Q_M)  % number of the wells
                r = dist_vec(i);
                if r<rw
                    r=rw;
                end

                R_influence = sqrt(2.246*perm*time/(visc_w*compr));        % pressure propagation radius for the time of injection [m]
                p_c = (Q_M(i,1)*0.159/365*12/86400*visc_w)/(2*pi*thick*perm)/1e6;        % characteristic pressure [MPa]

                % Baseline pressure increase by Q0
                Delta_p = Analytical_solution(r,R_influence,rc)*p_c;

                if time_count > 1
                    for v= 2:time_count
                        p_c = ((Q_M(i,v) - Q_M(i,v-1))*0.159/365*12/86400*visc_w)/(2*pi*thick*perm)/1e6;         % characteristic pressure [MPa]
                        R_influence = sqrt(2.246*perm*(time_step - (v-1)*dt)*86400*365/(visc_w*compr));       % pressure propagation radius for the time of injection [m]

                        Delta_p = Delta_p + Analytical_solution(r,R_influence,rc)*p_c;    % overpressure according to Nordbotten and Celia solution [MPa]
                    end
                end
                p_sup = p_sup +  Delta_p ;                             % superposed overpressure for the center wellbore [MPa]
            end

            p_fault_reservoir(fault_number,time_count) = p_sup ;        % pressure buildup on each fault for particular distance and well configuration
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
                for i = 1:length(Q_M)
                    r = dist_vec(i);
                    if r<rw
                        r=rw;
                    end

                    R_influence = sqrt(2.246*perm*time/(visc_w*compr));        % pressure propagation radius for the time of injection [m]
                    p_c = (Q_M(i,1)*0.159/365*12/86400*visc_w)/(2*pi*thick*perm)/1e6;        % characteristic pressure [MPa]

                    % Baseline pressure increase by Q0
                    Delta_p = Analytical_solution(r,R_influence,rc)*p_c;

                    if time_count > 1
                        for v= 2:time_count
                            p_c = ((Q_M(i,v) - Q_M(i,v-1))*0.159/365*12/86400*visc_w)/(2*pi*thick*perm)/1e6;         % characteristic pressure [MPa]
                            R_influence = sqrt(2.246*perm*(time_step - (v-1)*dt)*86400*365/(visc_w*compr));       % pressure propagation radius for the time of injection [m]

                            Delta_p = Delta_p + Analytical_solution(r,R_influence,rc)*p_c;    % overpressure [MPa]
                        end
                    end
                    p_sup = p_sup +  Delta_p ;
                end

                p_2Dgrid{time_count}(x,y) = p_sup ;               % pressure buildup across the reservoir for different times
            end
        end
    end
end
