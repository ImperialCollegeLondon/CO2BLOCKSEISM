function [deltap_cr,deltap_cr_ref,r_p,r_Sv,r_Shmin,r_SHmax,r_SHmax_dir,r_fault_mu,...
    r_fault_dip,r_fault_azi,pd_p,pd_Sv,pd_Shmin,pd_SHmax,pd_SHmax_dir,pd_fault_mu,...
    pd_fault_dip,pd_fault_azi] = Geomech_prob(fpath,fname,fault_dip,fault_azi,...
    f_dist_p,f_dist_Sv,f_dist_Shmin,f_dist_SHmax,f_dist_dir,f_dist_dip,f_dist_azi,...
    f_dist_mu,a_p_grad,a_Sv_grad,a_Shmin_grad,a_SHmax_grad,a_SHmax_dir,a_fault_dip,...
    a_fault_azi,a_fault_mu,nr_fault,n_MC)

    % -- read data
    [~,~,~,~,~,~,~,~,~,~,~,~,depth_seismic,depth_water,pres_grad_seismic,...
     Sv_grad_seismic,Shmin_grad_seismic,SHmax_grad_seismic,SHmax_dir_seismic,...
     fault_friction,pres_offset,Sv_offset,Shmin_offset,SHmax_offset] = read_data(fpath,fname);
    % --
        
    depth_seismic = depth_seismic - depth_water;

    % Distribution functions for uncertain variables
    if f_dist_p=='Uni'
        pd_p = makedist('Uniform','lower',pres_offset+depth_seismic*(pres_grad_seismic - a_p_grad/1000),...
            'upper',pres_offset+depth_seismic*(pres_grad_seismic + a_p_grad/1000)) ; 
    elseif f_dist_p=='Nor'
        pd_p = makedist('Normal','mu',pres_offset+depth_seismic*(pres_grad_seismic),...
            'sigma',depth_seismic*(a_p_grad/1000)) ; 
    end

    if f_dist_Sv=='Uni'
        pd_Sv = makedist('Uniform','lower',Sv_offset+depth_seismic*(Sv_grad_seismic - a_Sv_grad/1000),...
            'upper',Sv_offset+depth_seismic*(Sv_grad_seismic + a_Sv_grad/1000)) ; 
    elseif f_dist_Sv=='Nor'
        pd_Sv = makedist('Normal','mu',Sv_offset+depth_seismic*(Sv_grad_seismic),...
            'sigma',depth_seismic*(a_Sv_grad/1000)) ; 
    end

    if f_dist_Shmin=='Uni'
        pd_Shmin = makedist('Uniform','lower',Shmin_offset+depth_seismic*(Shmin_grad_seismic - a_Shmin_grad/1000),...
            'upper',Shmin_offset+depth_seismic*(Shmin_grad_seismic + a_Shmin_grad/1000)) ; 
    elseif f_dist_Shmin=='Nor'
        pd_Shmin = makedist('Normal','mu',Shmin_offset+depth_seismic*(Shmin_grad_seismic),...
            'sigma',depth_seismic*(a_Shmin_grad/1000)) ; 
    end
     
    if f_dist_SHmax=='Uni'
        pd_SHmax = makedist('Uniform','lower',SHmax_offset+depth_seismic*(SHmax_grad_seismic - a_SHmax_grad/1000),...
            'upper',SHmax_offset+depth_seismic*(SHmax_grad_seismic + a_SHmax_grad/1000)) ;  
    elseif f_dist_SHmax=='Nor'
        pd_SHmax = makedist('Normal','mu',SHmax_offset+depth_seismic*(SHmax_grad_seismic),...
            'sigma',depth_seismic*(a_SHmax_grad/1000)) ; 
    end

    if f_dist_dir=='Uni'
        pd_SHmax_dir = makedist('Uniform','lower',(SHmax_dir_seismic - a_SHmax_dir),...
            'upper',(SHmax_dir_seismic + a_SHmax_dir)) ;   
    elseif f_dist_dir=='Nor'
        pd_SHmax_dir = makedist('Normal','mu',SHmax_dir_seismic,...
            'sigma',a_SHmax_dir) ; 
    end

    if f_dist_mu=='Uni'
        pd_fault_mu = makedist('Uniform','lower',(fault_friction - a_fault_mu),...
            'upper',(fault_friction + a_fault_mu)) ; 
    elseif f_dist_mu=='Nor'
        pd_fault_mu = makedist('Normal','mu',fault_friction,...
            'sigma',a_fault_mu) ; 
    end
     

    % Random distribution of uncertain input parameters
    r_p = random(pd_p , n_MC , 1) ;
    r_Sv = random(pd_Sv , n_MC , 1) ;
    r_Shmin = random(pd_Shmin , n_MC , 1) ;
    r_SHmax = random(pd_SHmax , n_MC , 1) ;
    r_SHmax_dir = random(pd_SHmax_dir , n_MC , 1) ;
    r_fault_mu = random(pd_fault_mu , n_MC , 1) ;

    % Distribution functions and random generation of fault dip/azimuth
    % attributes
    for j= 1:nr_fault
        if f_dist_dip=='Uni'
            pd_fault_dip{j} = makedist('Uniform','lower',(fault_dip(j) - a_fault_dip),...
                'upper',(fault_dip(j) + a_fault_dip)) ;
        elseif f_dist_dip=='Nor'
            pd_fault_dip{j} = makedist('Normal','mu',fault_dip(j),...
                'sigma',a_fault_dip) ; 
            % limiting values to the range [0 90]
            pd_fault_dip{j} = truncate(pd_fault_dip{j},0,90);
        end
 
        if f_dist_azi=='Uni'
            pd_fault_azi{j} = makedist('Uniform','lower',(fault_azi(j) - a_fault_azi),...
                'upper',(fault_azi(j) + a_fault_azi)) ;
        elseif f_dist_azi=='Nor'
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
    
    % Calculating critical pressure for different combination of uncertain
    % parameters (Monte Carlo simultion)
    for i= 1:n_MC
        % calculating stress distribution along each fault
        for j= 1:nr_fault
            % for the number of Monte Carlo realizations
            [Sigma_n,Tau] = stress_projection(r_SHmax_dir(i),...
            r_SHmax(i),r_Shmin(i),r_Sv(i),r_fault_azi{j}(i),r_fault_dip{j}(i));

            deltap_cr{j}(i) = Sigma_n - r_p(i) - Tau/r_fault_mu(i);   % critical pressure to cause slip for MC realizations
        end
    end

    % Calculating critical pressure for all faults at reference stress
    % conditions and fault properties
      
    for i= 1:nr_fault
        [Sigma_n_ref,Tau_ref] = stress_projection(SHmax_dir_seismic,...
            (SHmax_offset + depth_seismic*SHmax_grad_seismic),...
            (Shmin_offset + depth_seismic*Shmin_grad_seismic),...
            (Sv_offset + depth_seismic*Sv_grad_seismic),...
            fault_azi(i),fault_dip(i));
        
        deltap_cr_ref(i) = Sigma_n_ref - (pres_offset + depth_seismic*pres_grad_seismic) - Tau_ref/fault_friction;   % critical pressure to cause slip for MC realizations
    end

end
