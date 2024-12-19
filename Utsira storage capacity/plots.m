%% Input data
w_id= 9;                        % wellbore configuration
n_distance = 2;                 % number of the desired interwell distance in the d_list. 
                                % Note that the maximum interwell distance for each scenario may not necessarily reside in the list 
time = 50;                      % time [year]
                                % if the selected time is not in the time series 0:dt:time_project, the closest time will be selected

axis_title_size = 9;           % font size of the axis title
title_size = 9;                % font size of the figure title
colorbar_size = 9;             % font size of the colorbar
axis_size = 9;                 % font size of numbers on the axes

plot_color = [1 1 1];           % Background plot color

% different scientific colormaps from crameri (alternatives: viridis, Jet)
cmap_fault = viridis(250);
cmap_contour = crameri('lajolla');
cmap_well = crameri('lajolla');

well_locations= 'Y';            % Y: plot      N: do not plot

% Histogram plot parameters
fault_No = 685;                 % number of the fault for which histograms are plotted
nr_bins = 40;                   % number of bins for each histogram

% seismicity control on injection
slip_prob_th = 0.5;             % Slip threshold to be considered as a potential seismic event
Mag_th = 4;                     % magnitude threshold for seismicity risk assessment
delta_sigma = 1;                % stress drop following fault slip [MPa]
nr_bins_FM = 20;                % number of bins for frequency-magnitude plot

% section "Seismicity magnitude calculation (probabilistic)"
nr_real = 1000;                  % number of realizations for event distribution analyses in probabilistic sense

if d_list(n_distance) > d_max_square(w_id)
    disp(sprintf('The selected interwell distance is not feasible \nOnly histogram can be plotted'))
end

save_figures= 'Y';              % saving plotted figures Yes: 'Y' or No: 'N'

resolution = 300;

fpath = pwd;

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

Utsira_boarders = xlsread("Utsira boundary points.xlsx");       % Utsira formation boundary
Utsira_south_boarders = xlsread("South Utsira points.xlsx");    % South Utsira boundary
Fault_data = xlsread("Faults.xlsx");                  % Fault attributes

% Assigning fault attributes (they are known in this case)
fault_dip = Fault_data(:,5);
fault_azi = Fault_data(:,4);

% fault length [m]
fault_length = Fault_data(:,3);

% fault coordinate (centroid) [m] 
fault_coord_x = Fault_data(:,1) - ref_X;
fault_coord_y = Fault_data(:,2) - ref_Y;

nr_fault = length(fault_dip);

% -- read data
[~,area_res,~,~,~,~,~,~,~,~,~,~,depth_seismic,depth_water,pres_grad_seismic,...
 Sv_grad_seismic,Shmin_grad_seismic,SHmax_grad_seismic,SHmax_dir_seismic,...
 fault_friction,pres_offset,Sv_offset,Shmin_offset,SHmax_offset] = read_data(fpath,fname);
    
%% Contour plots of pore pressure changes in response to CO2 injection (each plot in Fig. S5)
% for the selected well scenario, interwell distance and time
if d_list(n_distance) <= d_max_square(w_id)
    deltap_contour = figure;
    deltap_contour.Color = plot_color;
    deltap_contour.Units = 'centimeters';
    deltap_contour.Position = [6 8 9 7.2];    % [x y w h]
    
    time_series= dt:dt:time_project;

    [~,time_step] = min(abs(time_series-time));
    
    surf((Mesh_grid{w_id}(:,:,1))/1000 , (Mesh_grid{w_id}(:,:,2))/1000,p_2Dgrid{w_id}{n_distance}{time_step},'EdgeColor','none')
    grid off
    
    cbar = colorbar('Location','Eastoutside',...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

    set(get(cbar,'ylabel'),'string',['$\Delta$','{\it p}',' [MPa]'],...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);

    colormap(cmap_contour)

    clim([0 max(max(p_2Dgrid{w_id}{n_distance}{time_step}))]);
    
    hold on
    
    fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
    fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
    fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
    fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

    for i = 1:nr_fault
        plot3(fault_points_x(i,:) , fault_points_y(i,:),20+0*fault_points_x(i,:),...
            'LineWidth',0.5,'Color',[0.5 0.5 0.5])
        
        hold on
    end
    
    x_grid_num = well_list(2,w_id);
    y_grid_num = well_list(3,w_id);

    if well_locations == 'Y'
        % showing the position of wellbores
        distance = d_list(n_distance) * 1000 ;

        wells_coord_x_contour = (model_width - sqrt(area_res))*1000/2 + repmat((0:distance:distance*x_grid_num-1) + (y_grid_num*1000*d_max_square(w_id)-(x_grid_num-1)*distance)/2,y_grid_num,1);
        wells_coord_y_contour = (model_height - sqrt(area_res))*1000/2 + repmat(transpose((distance*(y_grid_num-1):-distance:0) + (y_grid_num*1000*d_max_square(w_id)-(y_grid_num-1)*distance)/2),1,x_grid_num);
  
        scatter3(wells_coord_x_contour/1000 , wells_coord_y_contour/1000,40*ones(well_list(3,w_id),well_list(2,w_id)), 10, ...
            '^','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',0.5)
    end
    hold on

    plot3((Utsira_south_boarders(:,1)-ref_X)/1000 , (Utsira_south_boarders(:,2)-ref_Y)/1000,...
        20+0*Utsira_south_boarders(:,1),'LineWidth',1.5,'Color',[0 0 0])
        
    text(67, 10, 6,'South Utsira', 'Interpreter','latex','FontSize', 10)
    
    title(['{\it t}','= ',num2str(time_series(time_step)),' y'],'FontSize', title_size)
    
    ylabel('Y, Northting [km]', ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);
    xlabel('X, Easting [km]', ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);

    shading interp
    view(2)
    axis('equal')

    axis([0 model_width 0 model_height])

    disp('contour plots of injection overpressure are done')
end


%% Histograms of uncertain parameters (Fig. S2)
Histogram_plot = figure;
Histogram_plot.Color = plot_color;
Histogram_plot.Units = 'centimeters';
Histogram_plot.Position = [6 8 19 7];    % [x y w h]

% p0
subplot(2,5,1);
hist_p = histogram(r_p,nr_bins);

xlabel(['{\it p}', ' [MPa]'],'Interpreter', 'latex','FontSize', axis_size)
label_p=ylabel('Realizations', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_size);

yticks([0 50 100 150])

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% Shmin
subplot(2,5,2);
hist_Sh = histogram(r_Shmin,nr_bins);
xlabel('$\sigma_h$ [MPa]','Interpreter', 'latex','FontSize', axis_size)

yticks([0 250 500])

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% SHmax
subplot(2,5,3);
hist_SH = histogram(r_SHmax,nr_bins);
xlabel('$\sigma_H$ [MPa]','Interpreter', 'latex','FontSize', axis_size)

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% Sv
subplot(2,5,4);
hist_Sv = histogram(r_Sv,nr_bins);
xlabel('$\sigma_v$ [MPa]','Interpreter', 'latex','FontSize', axis_size)

% SHmax azimuth
subplot(2,5,5);
hist_SHdir = histogram(r_SHmax_dir,nr_bins);
xlabel(['$\sigma_H$',' azimuth ', '[','$^\circ$',']'],'Interpreter', 'latex','FontSize', axis_size)

yticks([0 250 500])
set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% mu
subplot(2,5,6);
hist_mu = histogram(r_fault_mu,nr_bins);

xlabel('$\mu$','Interpreter', 'latex','FontSize', axis_size)
label_mu=ylabel('Realizations', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_size);

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% Fault stike
subplot(2,5,7);
hist_azi = histogram(r_fault_azi{fault_No},nr_bins);
xlabel(['Fault strike No ',num2str(fault_No),' [','$^\circ$',']'],'Interpreter', 'latex','FontSize', axis_size)

yticks([0 250 500])
set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% Fault dip
subplot(2,5,8);
hist_dip = histogram(r_fault_dip{fault_No},nr_bins);
xlabel(['Fault dip No ',num2str(fault_No),' [','$^\circ$',']'],'Interpreter', 'latex','FontSize', axis_size)

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% critical p
subplot(2,5,[9,10]);
hist_pc = histogram(deltap_cr{fault_No},nr_bins);
xlabel(['$\Delta$','{\it p}',' to slip on fault No ',num2str(fault_No), ' [MPa]'],'Interpreter', 'latex','FontSize', axis_size)

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

disp('histograms of probabilistic analyses of fault stability is done')


%% Torando sensitivity plot (Fig. S8)
Tornado_plot = figure;
Tornado_plot.Color = plot_color;
Tornado_plot.Units = 'centimeters';
Tornado_plot.Position = [6 8 8.38 7];    % [x y w h]

var_names={['{\it p}'];
'$\sigma_h$';
'$\sigma_H$';
'$\sigma_v$';
['Azimuth of ','$\sigma_H$'];
'$\mu$';
'Fault strike';
'Fault dip';
''};

tornado_baseline = deltap_cr_ref(fault_No);

% assigning sensitivity range for each parameter
% first column showing the minimum value of the variable
% second column showing the maximum value of the variable
if f_dist_p=='Uni'
    sens_range(1,1) = pres_offset + (depth_seismic - depth_water)*(pres_grad_seismic - a_p_grad/1000);
    sens_range(1,2) = pres_offset + (depth_seismic - depth_water)*(pres_grad_seismic + a_p_grad/1000);
elseif f_dist_p=='Nor'
    sens_range(1,1) = pres_offset + (depth_seismic - depth_water)*(pres_grad_seismic - 3*a_p_grad/1000);
    sens_range(1,2) = pres_offset + (depth_seismic - depth_water)*(pres_grad_seismic + 3*a_p_grad/1000);
end

if f_dist_Shmin=='Uni'
    sens_range(2,1) = Shmin_offset + (depth_seismic - depth_water)*(Shmin_grad_seismic - a_Shmin_grad/1000);
    sens_range(2,2) = Shmin_offset + (depth_seismic - depth_water)*(Shmin_grad_seismic + a_Shmin_grad/1000);
elseif f_dist_Shmin=='Nor'
    sens_range(2,1) = Shmin_offset + (depth_seismic - depth_water)*(Shmin_grad_seismic - 3*a_Shmin_grad/1000);
    sens_range(2,2) = Shmin_offset + (depth_seismic - depth_water)*(Shmin_grad_seismic + 3*a_Shmin_grad/1000);
end

if f_dist_SHmax=='Uni'
    sens_range(3,1) = SHmax_offset + (depth_seismic - depth_water)*(SHmax_grad_seismic - a_SHmax_grad/1000);
    sens_range(3,2) = SHmax_offset + (depth_seismic - depth_water)*(SHmax_grad_seismic + a_SHmax_grad/1000);
elseif f_dist_SHmax=='Nor'
    sens_range(3,1) = SHmax_offset + (depth_seismic - depth_water)*(SHmax_grad_seismic - 3*a_SHmax_grad/1000);
    sens_range(3,2) = SHmax_offset + (depth_seismic - depth_water)*(SHmax_grad_seismic + 3*a_SHmax_grad/1000);
end

if f_dist_Sv=='Uni'
    sens_range(4,1) = Sv_offset + (depth_seismic - depth_water)*(Sv_grad_seismic - a_Sv_grad/1000);
    sens_range(4,2) = Sv_offset + (depth_seismic - depth_water)*(Sv_grad_seismic + a_Sv_grad/1000);
elseif f_dist_Sv=='Nor'
    sens_range(4,1) = Sv_offset + (depth_seismic - depth_water)*(Sv_grad_seismic - 3*a_Sv_grad/1000);
    sens_range(4,2) = Sv_offset + (depth_seismic - depth_water)*(Sv_grad_seismic + 3*a_Sv_grad/1000);
end

if f_dist_dir=='Uni'
    sens_range(5,1) = SHmax_dir_seismic - a_SHmax_dir;
    sens_range(5,2) = SHmax_dir_seismic + a_SHmax_dir;
elseif f_dist_dir=='Nor'
    sens_range(5,1) = SHmax_dir_seismic - 3*a_SHmax_dir;
    sens_range(5,2) = SHmax_dir_seismic + 3*a_SHmax_dir;
end

if f_dist_mu=='Uni'
    sens_range(6,1) = fault_friction - a_fault_mu;
    sens_range(6,2) = fault_friction + a_fault_mu;
elseif f_dist_mu=='Nor'
    sens_range(6,1) = fault_friction - 3*a_fault_mu;
    sens_range(6,2) = fault_friction + 3*a_fault_mu; 
end

if f_dist_azi=='Uni'
    sens_range(7,1) = fault_azi(fault_No) - a_fault_azi;
    sens_range(7,2) = fault_azi(fault_No) + a_fault_azi;
elseif f_dist_azi=='Nor'
    sens_range(7,1) = fault_azi(fault_No) - 3*a_fault_azi;
    sens_range(7,2) = fault_azi(fault_No) + 3*a_fault_azi; 
end

if f_dist_dip=='Uni'
    sens_range(8,1) = fault_dip(fault_No) - a_fault_dip;
    sens_range(8,2) = fault_dip(fault_No) + a_fault_dip;
    if sens_range(8,2)>90
        sens_range(8,2) = 90;
    end
elseif f_dist_dip=='Nor'
    sens_range(8,1) = fault_dip(fault_No) - 3*a_fault_dip;
    sens_range(8,2) = fault_dip(fault_No) + 3*a_fault_dip; 
    if sens_range(8,2)>90
        sens_range(8,2) = 90;
    end
end

for i=1:length(var_names)-1
    % average values of the variable (baseline values)
    % i=1    pore pressue
    % i=2    Sh min
    % i=3    SH max
    % i=4    Sv
    % i=5    SH direction
    % i=6    Fault friction
    % i=7    Fault azimuth
    % i=8    Fault dip
    sens = [pres_offset + (depth_seismic - depth_water)*pres_grad_seismic;...
        Shmin_offset + (depth_seismic - depth_water)*Shmin_grad_seismic;...
        SHmax_offset + (depth_seismic - depth_water)*SHmax_grad_seismic;...
        Sv_offset + (depth_seismic - depth_water)*Sv_grad_seismic;...
        SHmax_dir_seismic;...
        fault_friction;...
        fault_azi(fault_No);...
        fault_dip(fault_No)];

    % assigning minimum value for the variable i
    sens(i) = sens_range (i,1); 

    [Sigma_n_ref,Tau_ref] = stress_projection(sens(5),sens(3),sens(2),sens(4),sens(7),sens(8));

    sens_deltap_low(i) = Sigma_n_ref - sens(1) - Tau_ref/sens(6);   % critical pressure to cause slip for minimum value of variable i
    
    % assigning maximum value for the variable i
    sens(i) = sens_range (i,2); 

    [Sigma_n_ref,Tau_ref] = stress_projection(sens(5),...
    sens(3),sens(2),sens(4),sens(7),sens(8));

    sens_deltap_high(i) = Sigma_n_ref - sens(1) - Tau_ref/sens(6);   % critical pressure to cause slip for maximum value of variable i

end
    sens_deltap_low(length(var_names)) = tornado_baseline;
    sens_deltap_high(length(var_names)) = tornado_baseline;

    % full range of variation of delta_p for each parameter 
    sens_deltap = abs(sens_deltap_high - sens_deltap_low);

    % sorting the parameters from large to small
    [sens_deltap,idx] = sort(sens_deltap(1:length(var_names)-1),'ascend');
    sens_deltap_high = [sens_deltap_high(idx) sens_deltap_high(9)];
    sens_deltap_low = [sens_deltap_low(idx) sens_deltap_low(9)];
    var_names = [var_names(idx)];

    h_tornado_low = barh(sens_deltap_low,'LineWidth',0.2);
    hold on
    h_tornado_high = barh(sens_deltap_high,'LineWidth',0.2);

    set(get(h_tornado_low,'BaseLine'),'BaseValue',tornado_baseline);

    title(['Sensitivity analysis for fault No ',num2str(fault_No)],'Interpreter', 'latex','FontSize', title_size)
    yticklabels(var_names);
    set(gca,'TickLabelInterpreter','latex')
    xlabel('Critical pressure for slip [MPa]', ...
            'Interpreter', 'latex','FontSize', axis_title_size);

    xlim([-2 10]);

    legend ([h_tornado_low h_tornado_high], {'Decrease of the input variable',...
        'Increase of the input variable'},...
        'Location','northwest','Color',[1 1 1],'box','on');

     set(gca,'Position',[0.265 0.13 0.7 0.81],'LineWidth',0.2,...
         'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

     disp('Tornado sensitivity plot is done')

%% Fault slip probability curves (Fig. 7a)
% only for the selected well scenario and interwell distance
if d_list(n_distance) <= d_max_square(w_id)
    cdf_plot = figure;
    cdf_plot.Color = plot_color;
    cdf_plot.Units = 'centimeters';
    cdf_plot.Position = [6 8 8.38 7];    % [x y w h]
    
    time_series= dt:dt:time_project;
    [~,time_step] = min(abs(time_series-time));
    
    for i=1:nr_fault
        [cdf x] = ecdf(deltap_cr{i});  
    
        % Matrix of the probability of slip on all faults at different times
        % "slip_prob_each_t" is only calculated for a certain "w_id" and "n_distance"
        for j = 1:length(time_series)
            [d , idx_prob] = min(abs(x - p_fault{w_id}{n_distance}(i,j)));
            slip_prob_each_t (i,j) = cdf(idx_prob);
    
            if p_fault{w_id}{n_distance}(i,j) < min(x)
                slip_prob_each_t (i,j) = 0;
            elseif p_fault{w_id}{n_distance}(i,j) > max(x)
                slip_prob_each_t (i,j) = 1;
            end
        end
    
        x_cdf(:,i) = x;                          % critical pressure arrays for all faults
        slip_prob_plot(:,i) = cdf;               % slip probability arrays for all faults
    end
    
    % index for assigning colors to each curve
    color_idx = linspace(0,1,length(cmap_fault));

    count = 0;
    order = [];
    for i=1:nr_fault
        if (fault_coord_x(i) >= -0.02*model_width*1000 ...
                && fault_coord_x(i) <= 1.02*model_width*1000) ...
                && (fault_coord_y(i) >= -0.02*model_height*1000 ...
                && fault_coord_y(i) <= 1.02*model_height*1000)
            
            count=count+1;
            
            [d , idx] = min(abs(color_idx - slip_prob_each_t (i,time_step)));
            order(count,1) = i;
            order(count,2) = idx;

        end
    end

    % plotting the probability curves in a way the most likely to slip
    % locate on top of the others
    order_sorted = sortrows(order,2);
    for i=1:count
        plot(x_cdf(:,order_sorted(i,1)),slip_prob_plot(:,order_sorted(i,1)),'LineWidth',1,'color',cmap_fault(order_sorted(i,2),:));
        hold on
    end
    
    title('')
    ylabel('Fault slip probability', ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);
    xlabel('Pore pressure changes [MPa]', ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);
    box('on');

    cbar = colorbar('Position',[0.875 0.13 0.0381 0.8], ...
        'Ticks',[0 0.2 0.4 0.6 0.8 1],... 
        'LineWidth',0.2,...
        'FontName', 'Arial', ...
        'TickLabelInterpreter','latex',...
        'FontSize',colorbar_size);
            
    set(get(cbar,'ylabel'),'string',['Fault slip probability (',...
    '{\it t}','= ',num2str(time_series(time_step)),' y)'],...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);
    caxis([0 1]);
    colormap(cmap_fault)

    xlim([0 10]);
    
    set(gca,'Position',[0.11 0.13 0.7 0.8],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

    hold off

    disp('fault slip probability curves are derived and plotted')
end


%% Spatial distribution of the faults colorcoded with slip probability at time t (Fig. 7b)
% only for the selected well scenario and interwell distance

if d_list(n_distance) <= d_max_square(w_id)   
    
    fault_slip_prob_plot = figure;
    fault_slip_prob_plot.Color = plot_color;
    fault_slip_prob_plot.Units = 'centimeters';
    fault_slip_prob_plot.Position = [6 8 8.38 7];    % [x y w h]
    
    axes('LineWidth',0.2);

    x_grid_num = well_list(2,w_id);
    y_grid_num = well_list(3,w_id);

    time_series= dt:dt:time_project;
    [~,time_step] = min(abs(time_series-time));
    
    % Fault positions colorcoded with the slip probability at time t
    
    color_idx = linspace(0,1,length(cmap_fault));                % index for assigning colors to each curve
    
    fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
    fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
    fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
    fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

    for i = 1:nr_fault
        [d , idx] = min(abs(color_idx -slip_prob_each_t(i,time_step)));
        
        plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',1.5,'Color',cmap_fault(idx,:))
    
        hold on
    end
    hold on

    plot((Utsira_south_boarders(:,1)-ref_X)/1000 , (Utsira_south_boarders(:,2)-ref_Y)/1000,...
        'LineWidth',1,'Color',[0 0 0])
    hold on
    
    if well_locations == 'Y'
        % wellbore locations
        distance = d_list(n_distance) * 1000 ;
        
        wells_coord_x = (model_width - sqrt(area_res))*1000/2 + repmat((0:distance:distance*x_grid_num-1) + (y_grid_num*1000*d_max_square(w_id)-(x_grid_num-1)*distance)/2,y_grid_num,1);
        wells_coord_y = (model_height - sqrt(area_res))*1000/2 + repmat(transpose((distance*(y_grid_num-1):-distance:0) + (y_grid_num*1000*d_max_square(w_id)-(y_grid_num-1)*distance)/2),1,x_grid_num);
                  
        scatter(wells_coord_x/1000 , wells_coord_y/1000, 7,'^','filled',...
            'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',0.5); hold on
    end
       
    text(67, 10,'South Utsira', 'Interpreter','latex','FontSize', 9)

    hold off
       
    box('on');
    
    cbar = colorbar('Position',[0.86 0.13 0.0381 0.8], ...
    'LineWidth',0.2,...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);
    
    set(get(cbar,'ylabel'),'string','Fault slip probability',...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);

    axis('equal')
        
    clim([0 1]);
    
    colormap(cmap_fault)
    
    title(['{\it t}', '= ',num2str(time), ' [y]'])
    ylabel('Y, Northting [km]', ...
        'Interpreter', 'latex', ...
        'FontSize', axis_title_size);
    xlabel('X, Easting [km]', ...
        'Interpreter', 'latex', ...
        'FontSize', axis_title_size);

    axis([0 model_width 0 model_height])

    xticks([0 30 60 90 120])
    yticks([0 30 60 90 120])

    set(gca,'Position',[0.06 0.13 0.8 0.8],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

    disp('distribution of fault slip probability at the specified time is plotted')
end  



%% Temporal evolution of fault slip probability (Fig. 10)
% only for the selected well scenario and interwell distance

if d_list(n_distance) <= d_max_square(w_id)   
    slipprob_temporal_fault = figure;
    slipprob_temporal_fault.Color = plot_color;
    slipprob_temporal_fault.Units = 'centimeters';
    slipprob_temporal_fault.Position = [6 8 8.38 7];    % [x y w h]
    
    time_series_plot= [0 dt:dt:time_project];
    
    count=0;
    for i = 1:nr_fault
        if (fault_coord_x(i) >= -0.02*model_width*1000 ...
            && fault_coord_x(i) <= 1.02*model_width*1000) ...
            && (fault_coord_y(i) >= -0.02*model_height*1000 ...
            && fault_coord_y(i) <= 1.02*model_height*1000)
            count=count+1;
            plot(time_series_plot,[0 slip_prob_each_t(i,:)],'LineWidth',1,...
                'Color',cmap_fault(order(count,2),:));
            hold on
        end
    end
    hold off
    
    cbar = colorbar('Position',[0.88 0.13 0.0381 0.85], ...
    'LineWidth',0.2,...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);
    
    set(get(cbar,'ylabel'),'string',['Fault slip probability (',...
        '{\it t}','= ',num2str(time_series(time_step)),' y)'],...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);
       
    clim([0 1]);
    
    colormap(cmap_fault)

    ylabel('Fault slip probability', ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);
    xlabel(['{\it t}',' [y]'], ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);

    set(gca,'Position',[0.11 0.13 0.7 0.85],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

    disp('temporal evolution of fault slip probability is plotted')
end



%% Seismicity magnitude calculation (probabilistic with threshold) 
% These plots account for seismic events until the selected time
% only for the selected well scenario and interwell distance
if d_list(n_distance) <= d_max_square(w_id)
    time_series= dt:dt:time_project;
    
    [~,time_step] = min(abs(time_series-time));

    for f=1:nr_fault
    %%%%%%%%%% Probabilistic with threshold
    % slip occurs if slip probability exceeds a user-assigned threshold
        if  max(slip_prob_each_t(f,1:time_step)) < slip_prob_th
            t_slip_prob(f) = NaN;
            Mag_prob_th(f) = NaN;
            cat_prob_coord_x(f) = NaN;         % predicted catalog from probability analysis with threshold
            cat_prob_coord_y(f) = NaN;         % predicted catalog from probability analysis with threshold

        else
            % finding the first time step at which slip_prob > slip_prob_threshold
            % (probabilistic approach)
            [~ , idx_prob] = find(((slip_prob_each_t(f,1:time_step)-slip_prob_th)>=0)==1);

            % time at which slip occurs on fault f (probabilistic approach) [y]
            % linear interpolation is considered between to consecutive pressure data 
            if min(idx_prob)>1
                t_slip_prob(f)= ((slip_prob_th - slip_prob_each_t(f,min(idx_prob)))/...
                    ((slip_prob_each_t(f,min(idx_prob)) - ...
                    slip_prob_each_t(f,min(idx_prob)-1))/(time_series(min(idx_prob))...
                    -time_series(min(idx_prob)-1))) + time_series(min(idx_prob)));
            else
                t_slip_prob(f)= 0.1;
            end

            % deltap at the slip time estimated in a probabilistic way
            if min(idx_prob)>1
                deltap_cr_prob(f) = ((p_fault{w_id}{n_distance}(f,min(idx_prob)) - ...
                    p_fault{w_id}{n_distance}(f,min(idx_prob)-1))/(time_series(min(idx_prob))...
                    -time_series(min(idx_prob)-1))) * (t_slip_prob(f)-time_series(min(idx_prob)))...
                    + p_fault{w_id}{n_distance}(f,min(idx_prob));
            else
                deltap_cr_prob(f) = p_fault{w_id}{n_distance}(f,1);
            end
      
            % largest possible earthquake magnitude
            % Circular fault model for moment calculation
            M0 = (16/7)*delta_sigma*1000000*(fault_length(f)/2)^3;      % moment [N.m]

            Mag_prob_th(f) = (2/3)*log10(M0) - 6.07;     % moment magnitude Mw 
                    
            cat_prob_coord_x(f) = fault_coord_x(f);
            cat_prob_coord_y(f) = fault_coord_y(f);
        end
    end
end


%% synthetic seismicity catalog (probabilistic with threshold analysis) (Fig. 9)
% only for the selected well scenario and interwell distance
if d_list(n_distance) <= d_max_square(w_id)

    catalog_syn = figure;
    catalog_syn.Color = plot_color;
    catalog_syn.Units = 'centimeters';
    catalog_syn.Position = [6 8 8.38 8.5];    % [x y w h]
    
    time_series= dt:dt:time_project;
    
    [~,time_step] = min(abs(time_series-time));
    
    ax1 = axes(catalog_syn,'LineWidth',0.2); 
        
    surfc((Mesh_grid{w_id}(:,:,1))/1000 , (Mesh_grid{w_id}(:,:,2))/1000,p_2Dgrid{w_id}{n_distance}{time_step},'EdgeColor','none')

    grid off

    box('on');
    
    cbar = colorbar('Position',[0.89 0.26 0.0381 0.7], ...
    'Ticks',[0 1 2 3 4],...
    'LineWidth',0.2,...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);
    
    set(get(cbar,'ylabel'),'string',['$\Delta$','{\it p}',' [MPa]'],...
    'Rotation',90,...
    'FontName', 'Arial', ...
    'interpreter', 'latex',...
    'fontsize',colorbar_size);
    
    colormap(cmap_contour)
    
    hold on

    fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
    fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
    fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
    fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

    for i = 1:nr_fault
        plot3(fault_points_x(i,:) , fault_points_y(i,:),20+0*fault_points_x(i,:),...
            'LineWidth',0.5,'Color',[0.5 0.5 0.5])
        
        hold on
    end

    x_grid_num = well_list(2,w_id);
    y_grid_num = well_list(3,w_id);

    if well_locations == 'Y'
        % showing the position of wellbores
        distance = d_list(n_distance) * 1000 ;

        wells_coord_x_contour = (model_width - sqrt(area_res))*1000/2 + repmat((0:distance:distance*x_grid_num-1) + (y_grid_num*1000*d_max_square(w_id)-(x_grid_num-1)*distance)/2,y_grid_num,1);
        wells_coord_y_contour = (model_height - sqrt(area_res))*1000/2 + repmat(transpose((distance*(y_grid_num-1):-distance:0) + (y_grid_num*1000*d_max_square(w_id)-(y_grid_num-1)*distance)/2),1,x_grid_num);
  
        scatter3(wells_coord_x_contour/1000 , wells_coord_y_contour/1000,40*ones(well_list(3,w_id),well_list(2,w_id)), 7, ...
            '^','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',0.5)
    end
    hold on

    plot3((Utsira_south_boarders(:,1)-ref_X)/1000 , (Utsira_south_boarders(:,2)-ref_Y)/1000,...
        20+0*Utsira_south_boarders(:,1),'LineWidth',1,'Color',[0 0 0])
        
    text(67, 10, 6,'South Utsira', 'Interpreter','latex','FontSize', 9)

    shading interp
    view(2)
    axis('equal')
    
    axis([0 model_width 0 model_height])

    xticks([0 30 60 90 120])
    yticks([0 30 60 90 120])
    
    title(['{\it t}','= ',num2str(time_series(time_step)),' y'],...
        'FontSize', title_size,'interpreter', 'latex')
    
    ylabel('Y, Northting [km]', ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);
    xlabel('X, Easting [km]', ...
                'Interpreter', 'latex', ...
                'FontSize', axis_title_size);

    ax2 = axes;
    
    % plotting seismicity data
    M3=scatter3(cat_prob_coord_x/1000 , cat_prob_coord_y/1000,15+0*cat_prob_coord_x,15,...
        Mag_prob_th,'o','filled','MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',0.5);
    
    clim([3.6 4.6]);

    cbar_mag = colorbar('Location','Southoutside',...
    'Ticks',[3.6 3.8 4 4.2 4.4 4.6],...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

    set(cbar_mag,'Position',[0.131 0.07 0.7 0.0381])
  
    set(get(cbar_mag,'ylabel'),'string','Seismic magnitude',...
    'FontName', 'Arial', ...
    'interpreter', 'latex',...
    'fontsize',colorbar_size);
    
    colormap(ax2,cmap_fault)
    
    ax2.UserData = linkprop([ax1,ax2],...
        {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
        'ydir','xdir','xlim','ylim'});
    
    shading interp
    view(2)
    axis('equal')
    
    axis([0 model_width 0 model_height])
    
    ax2.Visible = 'off';

    set(gca,'Position',[0.125 0.26 0.70 0.70],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

    text('String','a','Interpreter','latex','FontSize', 10,...
    'Position',[-20 123 0]);

    disp('plot for synthetic seismicity catalog is done')
end



%% Resources constrained by the maximum probability for M>threshold (Fig. 11)
time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

% matrix of maximum slip probability for at least one event > Mth (threshold)
% for all well scenarios and interwell distances at the specified time
scenarios_Mth_maxprob = zeros(width(well_list),nr_dist);

% matrix of the number of events (recognized by probability > threshold)
% larger than a threshold magnitude for all well scenarios and interwell 
% distances at the specified time
scenarios_Mth_number = zeros(width(well_list),nr_dist);

% matrix of maximum injection overpressure
% for all well scenarios and interwell distances at the specified time
scenarios_deltap = zeros(width(well_list),nr_dist);

% possible injection volume given the injection rate at the specified time
V_poss = zeros(width(well_list),nr_dist);

% calculating slip probability for each fault at the specified time
for i=1:width(well_list)         % all allowable well scenarios
    for j = 1:nr_dist
        % slip probability calculation only for allowable scenarios and distances
        if d_list(j) <= d_max_square(i)        % all allowable distances
            count_Mth_no = 0;
            prob_max = 0;
            count = 0;
            slip_probability = [];

            for f=1:nr_fault
                % largest possible earthquake magnitude
                % Circular fault model for moment calculation
                M0 = (16/7)*delta_sigma*1000000*(fault_length(f)/2)^3;      % moment [N.m]
                Mag = (2/3)*log10(M0) - 6.07;     % moment magnitude Mw 
                 
                if Mag >= Mag_th
                    count = count + 1;

                    [cdf x] = ecdf(deltap_cr{f});
                    [d , idx_prob] = min(abs(x - p_fault{i}{j}(f,time_step)));
    
                    slip_probability(count) = cdf(idx_prob);
    
                    if p_fault{i}{j}(f,time_step) < min(x)
                        slip_probability(count) = 0;
                    elseif p_fault{i}{j}(f,time_step) > max(x)
                        slip_probability(count) = 1;
                    end

                    % calculating maximum probability of at least one event
                    % with magnitude larger than a threshold
                    if slip_probability(count) > prob_max
                        % slip probability at the specified time
                        prob_max = slip_probability(count);
                    end

                    % calculating the number of events with magnitude
                    % larger than a threshold
                    % threshold approach for event recognition
                    if slip_probability(count) >= slip_prob_th
                        count_Mth_no = count_Mth_no + 1;
                    end
                end
            end


            scenarios_Mth_maxprob (i,j) = prob_max;
            scenarios_Mth_number (i,j) = count_Mth_no;
            scenarios_deltap (i,j) = max(max(p_2Dgrid{i}{j}{time_step}));
            % scenarios_Mth_number_probabilistic (i,j) = round(mean(sum(slip_realizations_Mth)));


%%%%%%%%%%% Calculations of injectable CO2 mass in each scenario %%%%%%%%%%
            % finding where t in the Q_M timing list is located
            bigger_value_location = find((Q_M(1,:)-time)>=0);
            if bigger_value_location>0
                previous_rate_no = bigger_value_location(1)-1;
            else
                [~,previous_rate_no] = size(Q_M);
            end

            % CO2 mass at a fixed injection rate over the specified time
            V_poss(i,j) = 0.001 * well_list(1,i) * ...
                ((Q_M(1,2:previous_rate_no)-Q_M(1,1:(previous_rate_no-1)))*...
                transpose(Q_M(2,1:previous_rate_no-1)) + ...
                (time - Q_M(1,previous_rate_no)) * Q_M(2,previous_rate_no));
        else
            % scenarios not complying with square reservoir shape
            scenarios_Mth_maxprob (i,j) = NaN;
            scenarios_Mth_number (i,j) = NaN;
            scenarios_deltap (i,j) = NaN; 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%% possible injection volume %%%%%%%%%%%%%%%%%%%%%%%%
% storage resources constrained by the maximum probability of at least one
% event with magnitude larger than a threshold
Resource_Mth_prob_plot = figure;
Resource_Mth_prob_plot.Color = plot_color;
Resource_Mth_prob_plot.Units = 'centimeters';
Resource_Mth_prob_plot.Position = [6 8 8.3 10.5];    % [x y w h]

d_mat = repmat(d_list,length(well_list),1);
d_max_square_mat = repmat(d_max_square',1,nr_dist);

% possible scenarios based on fault reactivation constraint (1:possible)
possible_faultstability = d_mat <= d_max_square_mat;

resource_Mth_prob = imagesc(scenarios_Mth_maxprob,'AlphaData',0.7*possible_faultstability);                 % colorcoded by the fault slip probability

cbar_Resource_Mth = colorbar('Location','Southoutside',...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

set(cbar_Resource_Mth,'Position',[0.12 0.07 0.86 0.0381])

set(get(cbar_Resource_Mth,'ylabel'),'string',...
    ['Maximum slip probability of M',' $\geq$ ', num2str(Mag_th)],...
    'FontName', 'Arial', ...
    'interpreter', 'latex',...
    'fontsize',colorbar_size);
    
colormap(cmap_fault)

clim([0 1]);

hold on;
% adding grid lines
for i = 1:length(well_list)
   plot([0.5,nr_dist+0.5],[i-0.5,i-0.5],...
       'LineWidth',0.1,'Color',[0 0 0]);
   Y_label{i} = sprintf('%.0f',well_list(1,i));
end

for i = 1:nr_dist
   plot([i-0.5,i-0.5],[0.5,length(well_list)+0.5],...
       'LineWidth',0.1,'Color',[0 0 0]);
   X_label{i} = sprintf('%.1f',d_list(i));
end

% Illustrating possible CO2 injection volume for each scenario on the cells
% based on constant injection rate over the specified period
for ii = 1:nr_dist
    for jj = 1:length(well_list)
        if V_poss(jj,ii)~=0
            if V_poss(jj,ii) < 10
                text(ii-0.3, jj+0.13, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            else if V_poss(jj,ii)>=10 && V_poss(jj,ii)<100
                text(ii-0.38, jj+0.11, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            else
                text(ii-0.47, jj+0.13, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            end
            end
        end
    end
end
hold off

% x and y axes labeled with the interwell distance and well numbers,
% respectively. Hide them to label by numbers
resource_Mth_prob.Parent.XTick = 1:nr_dist;
resource_Mth_prob.Parent.YTick = 1:length(well_list);

resource_Mth_prob.Parent.XTickLabel = X_label;
resource_Mth_prob.Parent.YTickLabel = Y_label;

title(['Injected CO','$_2$ ','volume [Gt]'])
    
ylabel('Number of the wells', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
xlabel('Interwell distance cases', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);

set(gca,'Position',[0.12 0.23 0.86 0.73],'LineWidth',0.2,'FontSize',axis_size,...
'TickDir','none','LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

text('String','a','Interpreter','latex','FontSize', 10,...
    'Position',[-1.3853550295858 -0.45580110497237 0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% storage resources constrained by the number of events larger than a
% threshold (being an event is recognized by a probability threshold)
Resource_Mth_NO_plot = figure;
Resource_Mth_NO_plot.Color = plot_color;
Resource_Mth_NO_plot.Units = 'centimeters';
Resource_Mth_NO_plot.Position = [6 8 8.3 10.5];    % [x y w h]

resource_Mth_NO = imagesc(scenarios_Mth_number,'AlphaData',0.7*possible_faultstability);                 % colorcoded by the fault slip probability

cbar_Resource_Mth_NO = colorbar('Location','Southoutside',...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

set(cbar_Resource_Mth_NO,'Position',[0.12 0.07 0.86 0.0381])

set(get(cbar_Resource_Mth_NO,'ylabel'),'string',...
    ['Number of events with M',' $\geq$ ', num2str(Mag_th)],...
    'FontName', 'Arial', ...
    'interpreter', 'latex',...
    'fontsize',colorbar_size);
    
colormap(cmap_fault)

hold on;
% adding grid lines
for i = 1:length(well_list)
   plot([0.5,nr_dist+0.5],[i-0.5,i-0.5],...
       'LineWidth',0.1,'Color',[0 0 0]);
   Y_label{i} = sprintf('%.0f',well_list(1,i));
end

for i = 1:nr_dist
   plot([i-0.5,i-0.5],[0.5,length(well_list)+0.5],...
       'LineWidth',0.1,'Color',[0 0 0]);
   X_label{i} = sprintf('%.1f',d_list(i));
end

% Illustrating possible CO2 injection volume for each scenario on the cells
% based on constant injection rate over the specified period
for ii = 1:nr_dist
    for jj = 1:length(well_list)
        if V_poss(jj,ii)~=0
            if V_poss(jj,ii) < 10
                text(ii-0.3, jj+0.13, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            else if V_poss(jj,ii)>=10 && V_poss(jj,ii)<100
                text(ii-0.38, jj+0.11, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            else
                text(ii-0.47, jj+0.13, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            end
            end
        end
    end
end

hold off

% x and y axes labeled with the interwell distance and well numbers,
% respectively. Hide them to label by numbers
resource_Mth_NO.Parent.XTick = 1:nr_dist;
resource_Mth_NO.Parent.YTick = 1:length(well_list);

resource_Mth_NO.Parent.XTickLabel = X_label;
resource_Mth_NO.Parent.YTickLabel = Y_label;

title(['Injected CO','$_2$ ','volume [Gt]'])
    
ylabel('Number of the wells', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
xlabel('Interwell distance cases', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);

set(gca,'Position',[0.12 0.23 0.86 0.73],'LineWidth',0.2,'FontSize',axis_size,...
'TickDir','none','LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

text('String','b','Interpreter','latex','FontSize', 10,...
    'Position',[-1.3853550295858 -0.45580110497237 0]);

disp('plots of storage capacity constrained by induced seismicity potential are done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% storage resources constrained by the injection overpressure (Fig. S6)
Resource_deltap_plot = figure;
Resource_deltap_plot.Color = plot_color;
Resource_deltap_plot.Units = 'centimeters';
Resource_deltap_plot.Position = [6 8 8.3 10.5];    % [x y w h]

resource_deltap = imagesc(scenarios_deltap,'AlphaData',0.7*possible_faultstability);                 % colorcoded by the fault slip probability

cbar_resource_deltap = colorbar('Location','Southoutside',...
    'Ticks',[0 1 2 3 4 5 6], ...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

clim([0 6])

set(cbar_resource_deltap,'Position',[0.12 0.07 0.86 0.0381])

set(get(cbar_resource_deltap,'ylabel'),'string',...
    ['$\Delta$','{\it p}',' [MPa]'],...
    'FontName', 'Arial', ...
    'interpreter', 'latex',...
    'fontsize',colorbar_size);
    
colormap(cmap_contour)

hold on;
% adding grid lines
for i = 1:length(well_list)
   plot([0.5,nr_dist+0.5],[i-0.5,i-0.5],...
       'LineWidth',0.1,'Color',[0 0 0]);
   Y_label{i} = sprintf('%.0f',well_list(1,i));
end

for i = 1:nr_dist
   plot([i-0.5,i-0.5],[0.5,length(well_list)+0.5],...
       'LineWidth',0.1,'Color',[0 0 0]);
   X_label{i} = sprintf('%.1f',d_list(i));
end

% Illustrating possible CO2 injection volume for each scenario on the cells
% based on constant injection rate over the specified period
for ii = 1:nr_dist
    for jj = 1:length(well_list)
        if V_poss(jj,ii)~=0
            if V_poss(jj,ii) < 10
                text(ii-0.3, jj+0.13, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            else if V_poss(jj,ii)>=10 && V_poss(jj,ii)<100
                text(ii-0.38, jj+0.11, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            else
                text(ii-0.47, jj+0.13, sprintf('%.1f',V_poss(jj,ii)), 'FontSize', 5.5);
            end
            end
        end
    end
end

hold off

% x and y axes labeled with the interwell distance and well numbers,
% respectively. Hide them to label by numbers
resource_deltap.Parent.XTick = 1:nr_dist;
resource_deltap.Parent.YTick = 1:length(well_list);

resource_deltap.Parent.XTickLabel = X_label;
resource_deltap.Parent.YTickLabel = Y_label;

title(['Injected CO','$_2$ ','volume [Gt]'])
    
ylabel('Number of the wells', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
xlabel('Interwell distance cases', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);

set(gca,'Position',[0.12 0.23 0.86 0.73],'LineWidth',0.2,'FontSize',axis_size,...
'TickDir','none','LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

disp('plot for maximum injection overpressure for different scenarios is done')

%% Saving figures
if save_figures=='Y'
    currentfolder = pwd;

    if isfolder([currentfolder,'\plots'])==0
        mkdir(currentfolder,'\plots');
    end

    filename=[currentfolder,'\plots\','deltap_contour'];
    print(deltap_contour, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(deltap_contour,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));

    filename=[currentfolder,'\plots\','Tornado_plot'];
    print(Tornado_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(Tornado_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));
    
    filename=[currentfolder,'\plots\','fault_slip_prob_plot'];
    print(fault_slip_prob_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(fault_slip_prob_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));

    filename=[currentfolder,'\plots\','cdf_plot'];
    print(cdf_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(cdf_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'))

    filename=[currentfolder,'\plots\','Histogram_plot'];
    print(Histogram_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(Histogram_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'))

    filename=[currentfolder,'\plots\','catalog_syn'];
    print(catalog_syn, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(catalog_syn,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'))

    filename=[currentfolder,'\plots\','Resource_Mth_prob_plot'];
    print(Resource_Mth_prob_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(Resource_Mth_prob_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'))

    filename=[currentfolder,'\plots\','Resource_Mth_NO_plot'];
    print(Resource_Mth_NO_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(Resource_Mth_NO_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'))

    filename=[currentfolder,'\plots\','Resource_deltap_plot'];
    print(Resource_deltap_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(Resource_deltap_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'))
end












