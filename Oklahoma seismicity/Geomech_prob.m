function [deltap_cr,deltap_cr_ref,r_p,r_Sv,r_Shmin,r_SHmax,r_SHmax_dir,r_fault_mu,...
    r_fault_dip,r_fault_azi,pd_p,pd_Sv,pd_Shmin,pd_SHmax,pd_SHmax_dir,pd_fault_mu,...
    pd_fault_dip,pd_fault_azi] = Geomech_prob(fpath,fname,stress_criticality,...
    dp_crit,fault_dip,fault_azi,f_dist_p,f_dist_Sv,f_dist_Shmin,f_dist_SHmax,...
    f_dist_dir,f_dist_dip,f_dist_azi,f_dist_mu,a_p_grad,a_Sv_grad,a_Shmin_grad,...
    a_SHmax_grad,a_SHmax_dir,a_fault_dip,a_fault_azi,a_fault_mu,nr_fault,n_MC)

    % -- read data
    [~,~,~,~,~,~,~,~,~,~,~,~,depth_seismic,pres_grad_seismic,...
     Sv_grad_seismic,Shmin_grad_seismic,SHmax_grad_seismic,SHmax_dir_seismic,...
     fault_friction] = read_data(fpath,fname);
    % --

    % Distribution functions for uncertain variables
    if strcmp(f_dist_p,'Uni')
        pd_p = makedist('Uniform','lower',depth_seismic*(pres_grad_seismic - a_p_grad/1000),...
            'upper',depth_seismic*(pres_grad_seismic + a_p_grad/1000)) ;
    elseif strcmp(f_dist_p,'Nor')
        pd_p = makedist('Normal','mu',depth_seismic*(pres_grad_seismic),...
            'sigma',depth_seismic*(a_p_grad/1000)) ;
    end

    if strcmp(f_dist_Sv,'Uni')
        pd_Sv = makedist('Uniform','lower',depth_seismic*(Sv_grad_seismic - a_Sv_grad/1000),...
            'upper',depth_seismic*(Sv_grad_seismic + a_Sv_grad/1000)) ;
    elseif strcmp(f_dist_Sv,'Nor')
        pd_Sv = makedist('Normal','mu',depth_seismic*(Sv_grad_seismic),...
            'sigma',depth_seismic*(a_Sv_grad/1000)) ;
    end

    if strcmp(f_dist_Shmin,'Uni')
        pd_Shmin = makedist('Uniform','lower',depth_seismic*(Shmin_grad_seismic - a_Shmin_grad/1000),...
            'upper',depth_seismic*(Shmin_grad_seismic + a_Shmin_grad/1000)) ; %#ok<LLMNC>
    elseif strcmp(f_dist_Shmin,'Nor')
        pd_Shmin = makedist('Normal','mu',depth_seismic*(Shmin_grad_seismic),...
            'sigma',depth_seismic*(a_Shmin_grad/1000)) ;
    end

    if stress_criticality == 'Y'
        pd_SHmax = [];
        q = ((1 + fault_friction^2)^0.5 + fault_friction)^2;
        % calculating reference SHmax in the case of stress criticality
        SHmax_grad_seismic = q*(Shmin_grad_seismic-pres_grad_seismic-dp_crit/depth_seismic) +...
            pres_grad_seismic + dp_crit/depth_seismic;
    else
        if strcmp(f_dist_SHmax,'Uni')
            pd_SHmax = makedist('Uniform','lower',depth_seismic*(SHmax_grad_seismic - a_SHmax_grad/1000),...
                'upper',depth_seismic*(SHmax_grad_seismic + a_SHmax_grad/1000)) ;  %#ok<LLMNC>
        elseif strcmp(f_dist_SHmax,'Nor')
            pd_SHmax = makedist('Normal','mu',depth_seismic*(SHmax_grad_seismic),...
                'sigma',depth_seismic*(a_SHmax_grad/1000)) ;
        end
    end

    if strcmp(f_dist_dir,'Uni')
        pd_SHmax_dir = makedist('Uniform','lower',(SHmax_dir_seismic - a_SHmax_dir),...
            'upper',(SHmax_dir_seismic + a_SHmax_dir)) ;
    elseif strcmp(f_dist_dir,'Nor')
        pd_SHmax_dir = makedist('Normal','mu',SHmax_dir_seismic,...
            'sigma',a_SHmax_dir) ;
    end

    if strcmp(f_dist_mu,'Uni')
        pd_fault_mu = makedist('Uniform','lower',(fault_friction - a_fault_mu),...
            'upper',(fault_friction + a_fault_mu)) ;
    elseif strcmp(f_dist_mu,'Nor')
        pd_fault_mu = makedist('Normal','mu',fault_friction,...
            'sigma',a_fault_mu) ;
    end


    % Random distribution of uncertain input parameters
    r_p = random(pd_p , n_MC , 1) ;
    r_Sv = random(pd_Sv , n_MC , 1) ;
    r_Shmin = random(pd_Shmin , n_MC , 1) ;
    r_SHmax_dir = random(pd_SHmax_dir , n_MC , 1) ;
    r_fault_mu = random(pd_fault_mu , n_MC , 1) ;
    if stress_criticality == 'Y'
        % Approach 1: using reference pressure and friction values
        % q = ((1 + fault_friction^2)^0.5 + fault_friction)^2;
        % r_SHmax = q*(r_Shmin-pres_grad_seismic*depth_seismic-dp_crit) +...
        %     pres_grad_seismic*depth_seismic + dp_crit;

        % Approach 2: only using reference friction value
        % q = ((1 + fault_friction^2)^0.5 + fault_friction)^2;
        % r_SHmax = q*(r_Shmin-r_p-dp_crit) + r_p + dp_crit;

        % Approach 3: all parameters are exactly those of realizations
        q = ((1 + r_fault_mu.^2).^0.5 + r_fault_mu).^2;
        r_SHmax = q.*(r_Shmin-r_p-dp_crit) + r_p + dp_crit;
    else
        r_SHmax = random(pd_SHmax , n_MC , 1) ;
    end

    pd_fault_dip = cell(1,nr_fault);
    pd_fault_azi = cell(1,nr_fault);
    r_fault_dip = cell(1,nr_fault);
    r_fault_azi = cell(1,nr_fault);
    % Distribution functions and random generation of fault dip/azimuth
    % attributes
    for j= 1:nr_fault
        if strcmp(f_dist_dip,'Uni')
            pd_fault_dip{j} = makedist('Uniform','lower',(fault_dip(j) - a_fault_dip),...
                'upper',(fault_dip(j) + a_fault_dip)) ;
        elseif strcmp(f_dist_dip,'Nor')
            pd_fault_dip{j} = makedist('Normal','mu',fault_dip(j),...
                'sigma',a_fault_dip) ;
            % limiting values to the range [0 90]
            pd_fault_dip{j} = truncate(pd_fault_dip{j},0,90);
        end

        if strcmp(f_dist_azi,'Uni')
            pd_fault_azi{j} = makedist('Uniform','lower',(fault_azi(j) - a_fault_azi),...
                'upper',(fault_azi(j) + a_fault_azi)) ;
        elseif strcmp(f_dist_azi,'Nor')
            pd_fault_azi{j} = makedist('Normal','mu',fault_azi(j),...
                'sigma',a_fault_azi) ;
        end

        % Matrices with nr_fault columns and n_MC rows
        r_fault_dip{j} = random(pd_fault_dip{j} , n_MC , 1) ;
        r_fault_azi{j} = random(pd_fault_azi{j} , n_MC , 1) ;

        % adjusting azimuth values >360 or <0
        for i= 1:n_MC
            if r_fault_azi{j}(i)>360
                r_fault_azi{j}(i) = r_fault_azi{j}(i) - 360;
            elseif r_fault_azi{j}(i)<0
                r_fault_azi{j}(i) = r_fault_azi{j}(i) + 360;
            end
        end
    end
    
    r_SHmax_dir_rad = deg2rad(r_SHmax_dir);
    r_fault_azi_rad = cellfun(@deg2rad,r_fault_azi,'UniformOutput',false);
    r_fault_dip_rad = cellfun(@deg2rad,r_fault_dip,'UniformOutput',false);
    deltap_cr = cell(1,nr_fault);
    % calculating stress distribution along each fault
    for j= 1:nr_fault
        % Calculating critical pressure for different combination of uncertain
        % parameters (Monte Carlo simultion)
        [Sigma_n,Tau] = arrayfun(@stress_projection,...
            r_SHmax_dir_rad,r_SHmax,r_Shmin,r_Sv,r_fault_azi_rad{j},r_fault_dip_rad{j});
        deltap_cr{j} = Sigma_n - r_p - Tau./r_fault_mu;
        if mod(j,100) == 0
            disp(j);
        end
    end

    deltap_cr_ref = zeros(1,nr_fault);
    % Calculating critical pressure for all faults at reference stress
    % conditions and fault properties
    for i= 1:nr_fault
        [Sigma_n_ref,Tau_ref] = stress_projection(SHmax_dir_seismic,...
            depth_seismic*SHmax_grad_seismic,depth_seismic*Shmin_grad_seismic,depth_seismic*Sv_grad_seismic,...
            fault_azi(i),fault_dip(i)); %#ok<LLMNC>
        % critical pressure to cause slip for MC realizations
        deltap_cr_ref(i) = Sigma_n_ref - depth_seismic*pres_grad_seismic ...
            - Tau_ref/fault_friction;
    end

end
