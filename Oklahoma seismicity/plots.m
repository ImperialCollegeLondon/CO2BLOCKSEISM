%% Input data

time = 18;                     % time [year] for plotting
                               % if the selected time is not in the time series 0:dt:time_project, the closest time will be selected

axis_title_size = 9;           % font size of the axis title
title_size = 9;                % font size of the figure title
colorbar_size = 9;             % font size of the colorbar
axis_size = 9;                 % font size of numbers on the axes

plot_color = [1 1 1];          % Background plot color

% Assigning colormaps
% different scientific colormaps from crameri can be used 
% Crameri, F. (2018). Scientific colour maps. Zenodo. at:
% https://doi.org/10.5281/zenodo.1243862

% (alternatives: viridis, Jet)

cmap_fault = viridis(256);          % colormap of fault slip probability

cmap_contour = crameri('lajolla');  % colormap of contour plot
cmap_well = crameri('lajolla');     % colormap of well distribution with their cumulative injected fluid volume

% assigning the range of the colorbar
max_color_pressure = 10;        % leave it empty [] for considering max possible value
max_color_probability = 0.3;    % leave it empty [] for considering max possible value

cdf_color = 'probability';      % either 'probability' or 'slip pressure'
                                % probability: probability of fault slip at time t
                                % slip pressure: pressure to cause slip for the reference stress and fault attributes

well_locations= 'Y';            % Y: plot      N: do not plot

% Histogram plot parameters
fault_No = 2;                   % number of the fault for which histograms are plotted 
                                % use 2 to reproduce Fig. S2
nr_bins = 40;                   % number of bins for each histogram

save_figures= 'N';              % saving plotted figures Yes: 'Y' or No: 'N'

resolution = 300;

fpath = pwd;

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% Oklahoma data
County_boarders = xlsread("OK_KS_boarders.xlsx");     % Coordinate of county boarders
OK_state_boarder = xlsread("OK_state_border.xlsx");   % Coordinate of the whole OK state border
Seismicity_M3 = xlsread("Seismicity M3.xlsx");        % Seismic events of M>3
Seismicity_M4_5 = xlsread("Seismicity M4.5.xlsx");    % Seismic events of M>4.5
M3_monthlyrate = xlsread("Seismicity M3-monthly rate.xlsx");    % Monthly rate of M>3
well_data = xlsread("Inj_rate.xlsx");                 % Monthly injection rate data of the wellbores 
Fault_data = xlsread("Faults.xlsx");                  % Fault attributes

% Oklahoma State boarder coordinates
OK_state_boarder_x = 111.320*1000*(OK_state_boarder(:,2) - ref_lon).*cos(deg2rad(OK_state_boarder(:,1)));            
OK_state_boarder_y = 110.574*1000*(OK_state_boarder(:,1) - ref_lat);

% Study area coordinates
Studyarea_x = 111.320*1000*(Study_area(:,2) - ref_lon).*cos(deg2rad(Study_area(:,1)));            
Studyarea_y = 110.574*1000*(Study_area(:,1) - ref_lat);

% Loading seismicity data
event_M3_coord_x = 111.320*1000*(Seismicity_M3(:,2) - ref_lon).*cos(deg2rad(Seismicity_M3(:,1)));            
event_M3_coord_y = 110.574*1000*(Seismicity_M3(:,1) - ref_lat);

event_M4_5_coord_x = 111.320*1000*(Seismicity_M4_5(:,2) - ref_lon).*cos(deg2rad(Seismicity_M4_5(:,1)));            
event_M4_5_coord_y = 110.574*1000*(Seismicity_M4_5(:,1) - ref_lat);


% Assigning fault attributes (they are known in this case)
fault_dip = Fault_data(:,5);
fault_azi = Fault_data(:,4);

% fault length [m]
fault_length = Fault_data(:,3);

% fault coordinate (centroid) [m] 
fault_coord_x = Fault_data(:,1);
fault_coord_y = Fault_data(:,2);

nr_fault = length(fault_dip);

% Loading injection rate data [bbl/month]
Inj_data = well_data(:,5:220);

% Cumulative injected volume [m3]
cum_water = well_data(:,4)*0.159;

%% Injection and seismicity data in Oklahoma (Fig. 2)
Monthly_inj_rate_plot = figure;
Monthly_inj_rate_plot.Color = plot_color;
Monthly_inj_rate_plot.Units = 'centimeters';
Monthly_inj_rate_plot.Position = [6 8 15 9];    % [x y w h]

time_series= 2000+(1/12:1/12:time_project);

Inj_rate_total = sum(Inj_data)*0.159;

Inj_rate_cum(1) = Inj_rate_total(1);
for i=1:length(Inj_rate_total)-1
    Inj_rate_cum(i+1)  = Inj_rate_cum(i) + Inj_rate_total(i+1);
end

axes1 = axes('Position',[0.18 0.13 0.72 0.8],'FontSize',axis_size);   % Innerposition [left,bottom,width,height]
hold(axes1,'on');
box(axes1,'on');

xlim([2000 2018]);

yyaxis left
plot(time_series,Inj_rate_total,'LineWidth',1)

ylabel('Injection rate [$m^3$/month]', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
xlabel('Year', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
ylim([0 2e7]);

yyaxis right
plot(M3_monthlyrate(:,1),M3_monthlyrate(:,2),'LineWidth',1)

ylabel(['Monthly number of events ' '{\it M}' ' $\geq$ 3'], ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
% Parague
annotation(Monthly_inj_rate_plot,'arrow',[0.656922222222222 0.656922222222221],...
    [0.524574074074074 0.459618055555557],'HeadWidth',6,'HeadLength',5);
% Fairview
annotation(Monthly_inj_rate_plot,'arrow',[0.781099999999999 0.816377777777774],...
    [0.883833333333334 0.897062500000004],'HeadWidth',6,'HeadLength',5);
% Pawnee
annotation(Monthly_inj_rate_plot,'arrow',[0.868588888888889 0.854477777777775],...
    [0.574270833333333 0.516062500000002],'HeadWidth',6,'HeadLength',5);
% Cushing
annotation(Monthly_inj_rate_plot,'arrow',[0.83331111111111 0.853066666666664],...
    [0.288520833333334 0.336145833333334],'HeadWidth',6,'HeadLength',5);

% Parague
text('String',['Parague',sprintf('\n'),'{\it M} = 5.7'],...
    'Position',[2011 64.24 0],'FontSize',8);
% Parague
text('String',['Fairview',sprintf('\n'),'{\it M} = 5.1'],...
    'Position',[2013.09 110.2 0],'FontSize',8);
% Pawnee
text('String',['Pawnee',sprintf('\n'),'{\it M} = 5.8'],...
    'Position',[2016.2 73.37 0],'FontSize',8);
% Cushing
text('String',['Pawnee',sprintf('\n'),'{\it M} = 5.0'],...
    'Position',[2015.02352941177 16.7326732673267 0],'FontSize',8);
hold(axes1,'off');

axes2 = axes('Position',[0.08 0.13 0.82 0.8],...    % Innerposition [left,bottom,width,height]
          'FontSize',axis_size,'Color','none','XColor','none', ...
          'YColor',[0.3,0.3,0.3],'XLim',[1997.4 2018],'YLim',[0 14e8],...
          'XTick',[],'XTickLabel',[]);

hold(axes2,'on');
box(axes2,'off');

plot(axes2,time_series,Inj_rate_cum,'LineWidth',1,'Color',[0.3,0.3,0.3])

ylabel('Cumulative injected volume [$m^3$]', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);

hold(axes2,'off');

% OK state boarder
axes3 = axes('Position',[0.2 0.4 0.4 0.6],...    % Innerposition [left,bottom,width,height]
          'FontSize',axis_size,'Color','none','XColor','none','YColor','none',...
          'XTick',[],'XTickLabel',[]);
ylim([-200 350]);

hold(axes3,'on');
box(axes3,'off');

plot(OK_state_boarder_x/1000 , OK_state_boarder_y/1000,'LineWidth',0.2,'Color',[0 0 0])
hold on

plot(Studyarea_x/1000 , Studyarea_y/1000,'LineWidth',0.2,'LineStyle','--','Color',[0 0 0])
hold on

% plotting seismicity data
scatter(event_M3_coord_x/1000 , event_M3_coord_y/1000,1,'o','filled','MarkerFaceColor',[0.85,0.33,0.10]);

text('String','Oklahoma','Position',[-193.18 191.75 0],'FontSize',8);
text('String','Kansas','Position',[-193.18 241.24 0],'FontSize',8);
text('String','Study area','Position',[65.27 301.73 0],'FontSize',8);

axis('equal')

disp('general injection and seismicity data in Oklahoma is plotted')


%% Wellbore distribution (Fig. 3a)
well_distribution = figure;
well_distribution.Color = plot_color;
well_distribution.Units = 'centimeters';
well_distribution.Position = [6 8 9 7.2];    % [x y w h]

axes('LineWidth',0.1);

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        'LineWidth',0.1,'Color',[0 0 0])
    hold on
end

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

% plotting the faults
for i = 1:nr_fault
    plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',0.5,'Color',[0.5 0.5 0.5])

    hold on
end

% plotting the wellbores
scatter(wells_coord_x/1000 , wells_coord_y/1000,7,cum_water,'^','filled')

axis('equal')

axis([0 model_width 0 model_height])

ylabel('Y, Northting [km]', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);

% cbar = colorbar('Position',[0.79 0.13 0.037 0.79],...

cbar = colorbar('Position',[0.81 0.129 0.0381 0.7967],... 
'LineWidth',0.2,...
'FontName', 'Arial', ...
'TickLabelInterpreter','latex',...
'FontSize',colorbar_size);

set(get(cbar,'ylabel'),'string','Cumulative injection volume [$m^{3}$]',...
'Rotation',90,...
'FontName', 'Arial', ...
'interpreter', 'latex',...
'fontsize',colorbar_size);

set(gca,'ColorScale','log')

box('on');
axis('equal')

text(216, 210,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,'Kansas', 'Interpreter','latex','FontSize', 8)

text(-50, 315,'a', 'Interpreter','latex','FontSize', 10)

colormap(cmap_well)
xlim([0 300]);
ylim([0 300]);

xticks([0 100 200 300])
yticks([0 100 200 300])

set(gca,'Position',[0.0506 0.13 0.775 0.795],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1);

disp('wellbore distribution plot is done')


%% Pressure contour plot (Fig. 3b)
deltap_contour = figure;
deltap_contour.Color = plot_color;
deltap_contour.Units = 'centimeters';
deltap_contour.Position = [6 8 9 7.2];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;

[~,time_step] = min(abs(time_series-time));

surfc((Mesh_grid(:,:,1))/1000 , (Mesh_grid(:,:,2))/1000,p_2Dgrid{time_step},'EdgeColor','none')
grid off

cbar = colorbar('Position',[0.81 0.129 0.0381 0.7967],... 
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

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot3(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        2+0*County_boarders(:,m),'LineWidth',0.1,'Color',[0 0 0])
    hold on
end

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

% plotting the faults
for i = 1:nr_fault
    plot3(fault_points_x(i,:) , fault_points_y(i,:),2+0*fault_points_x(i,:),...
        'LineWidth',0.5,'Color',[0.5 0.5 0.5])

    hold on
end


% plotting seismicity data
M3=scatter3(event_M3_coord_x/1000 , event_M3_coord_y/1000,15+0*event_M3_coord_x,7,'o','filled',...
    'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.5);
hold on
m45=scatter3(event_M4_5_coord_x/1000 , event_M4_5_coord_y/1000,25+0*event_M4_5_coord_x,...
    40,'pentagram','filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],'LineWidth',0.5);

legend ([M3 m45], {['{\it M}' ' $\geq$ 3'] ['{\it M}' ' $\geq$ 4.5']},...
    'Location','southwest','Color',[1 1 1],'box','on');  

title('Arbuckle Formation, Dec. 2017','FontSize', title_size,...
    'interpreter', 'latex')

ylabel('Y, Northting [km]', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
            'Interpreter', 'latex', ...
            'FontSize', axis_title_size);

box('on');
shading interp
view(2)
axis('equal')

text(216, 210,1,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,1,'Kansas', 'Interpreter','latex','FontSize', 8)
text(-50, 315,'b', 'Interpreter','latex','FontSize', 10)

axis([0 model_width 0 model_height])

xticks([0 100 200 300])
yticks([0 100 200 300])

set(gca,'Position',[0.0506 0.13 0.775 0.795],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

disp('contour plot of pressure buildup is done')


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

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% Shmin
subplot(2,5,2);
hist_Sh = histogram(r_Shmin,nr_bins);
xlabel('$\sigma_h$ [MPa]','Interpreter', 'latex','FontSize', axis_size)

% yticks([0 250 500])
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

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% SHmax azimuth
subplot(2,5,5);
hist_SHdir = histogram(r_SHmax_dir,nr_bins);
xlabel(['$\sigma_H$',' azimuth ', '[','$^\circ$',']'],'Interpreter', 'latex','FontSize', axis_size)

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

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% Fault dip
subplot(2,5,8);
hist_dip = histogram(r_fault_dip{fault_No},nr_bins);
xlabel(['Fault dip No',num2str(fault_No),' [','$^\circ$',']'],'Interpreter', 'latex','FontSize', axis_size)

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

% critical p
subplot(2,5,[9,10]);
hist_pc = histogram(deltap_cr{fault_No},nr_bins);
xlabel(['$\Delta$','{\it p}',' to slip on fault No ',num2str(fault_No), ' [MPa]'],'Interpreter', 'latex','FontSize', axis_size)

set(gca,'FontSize',axis_size,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

disp('histograms of uncertain geomechanical parameters are plotted')


%% Fault slip probability curves (Figs. 5 and S4)
cdf_plot = figure;
cdf_plot.Color = plot_color;
cdf_plot.Units = 'centimeters';
cdf_plot.Position = [6 8 10 8.3];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

for i=1:nr_fault
    [cdf x] = ecdf(deltap_cr{i});  

    % Matrix of the probability of slip on all faults at different times
    for j = 1:length(time_series)
        [d , idx_prob] = min(abs(x - p_fault(i,j)));
        slip_prob_each_t (i,j) = cdf(idx_prob);

        if p_fault(i,j) < min(x)
            slip_prob_each_t (i,j) = 0;
        elseif p_fault(i,j) > max(x)
            slip_prob_each_t (i,j) = 1;
        end
    end

    x_cdf(:,i) = x;                          % critical pressure arrays for all faults
    slip_prob_plot(:,i) = cdf;               % slip probability arrays for all faults
end

% index for assigning colors to each curve
if strcmp(cdf_color,'probability')
    if isempty (max_color_probability) == 1
        color_idx = linspace(0,1,length(cmap_fault));
    else
        color_idx = linspace(0,max_color_probability,length(cmap_fault));
    end
elseif strcmp(cdf_color,'slip pressure')
    if isempty (max_color_pressure) == 1
        color_idx = linspace(min(deltap_cr_ref),max(deltap_cr_ref),length(cmap_fault));
    else
        color_idx = linspace(min(deltap_cr_ref),max_color_pressure,length(cmap_fault));
    end
else
    disp('Select the cdf colorcode from the two available options');
end

for i=1:nr_fault
    if strcmp(cdf_color,'probability')
        [d , idx] = min(abs(color_idx - slip_prob_each_t (i,time_step)));
        plot(x_cdf(:,i),slip_prob_plot(:,i),'LineWidth',1,'color',cmap_fault(idx,:));
    elseif strcmp(cdf_color,'slip pressure')
        [d , idx] = min(abs(color_idx - deltap_cr_ref (i)));
        plot(x_cdf(:,i),slip_prob_plot(:,i),'LineWidth',1,'color',cmap_fault(idx,:));
    end
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

cbar = colorbar('Position',[0.86 0.11 0.0381 0.85], ...
    'Ticks',[0 0.05 0.1 0.15 0.2 0.25 0.3],... 
    'LineWidth',0.2,...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

if strcmp(cdf_color,'probability')
    set(get(cbar,'ylabel'),'string','Fault slip probability (Dec. 2017)',...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);

    if isempty (max_color_probability) == 1
        clim([0 1]);
    else
        clim([0 max_color_probability]);
    end
    colormap(cmap_fault)

elseif strcmp(cdf_color,'slip pressure')
    set(get(cbar,'ylabel'),'string',['$\Delta$','{\it p}',' to slip [MPa]'],...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);
    if isempty (max_color_pressure) == 1
        clim([min(deltap_cr_ref) max(deltap_cr_ref)]);
    else
        clim([min(deltap_cr_ref) max_color_pressure]);
    end
    colormap(cmap_fault)
end

xlim([0 120]);

set(gca,'Position',[0.1 0.11 0.7 0.85],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

hold off

cdf_plot_trunc = figure;
cdf_plot_trunc.Color = plot_color;
cdf_plot_trunc.Units = 'centimeters';
cdf_plot_trunc.Position = [6 8 10 8.3];    % [x y w h]

axes('LineWidth',0.2);

for i=1:nr_fault
    if strcmp(cdf_color,'probability')
        [d , idx] = min(abs(color_idx - slip_prob_each_t (i,time_step)));
        plot(x_cdf(:,i),slip_prob_plot(:,i),'LineWidth',1,'color',cmap_fault(idx,:));
    elseif strcmp(cdf_color,'slip pressure')
        [d , idx] = min(abs(color_idx - deltap_cr_ref (i)));
        plot(x_cdf(:,i),slip_prob_plot(:,i),'LineWidth',1,'color',cmap_fault(idx,:));
    end
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

cbar = colorbar('Position',[0.86 0.11 0.0381 0.85], ...
    'Ticks',[0 0.05 0.1 0.15 0.2 0.25 0.3],... 
    'LineWidth',0.2,...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

if strcmp(cdf_color,'probability')
    set(get(cbar,'ylabel'),'string','Fault slip probability (Dec. 2017)',...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);

    if isempty (max_color_probability) == 1
        clim([0 1]);
    else
        clim([0 max_color_probability]);
    end
    colormap(cmap_fault)

elseif strcmp(cdf_color,'slip pressure')
    set(get(cbar,'ylabel'),'string',['$\Delta$','{\it p}',' to slip [MPa]'],...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'latex',...
        'fontsize',colorbar_size);
    if isempty (max_color_pressure) == 1
        clim([min(deltap_cr_ref) max(deltap_cr_ref)]);
    else
        clim([min(deltap_cr_ref) max_color_pressure]);
    end
    colormap(cmap_fault)
end

xlim([0 1]);

set(gca,'Position',[0.1 0.11 0.7 0.85],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

hold off

disp('calculations and plots of fault slip probability is done')

%% Spatial distribution of the faults colorcoded with slip probability at time t (Fig. 6a)
fault_slip_prob_plot = figure;
fault_slip_prob_plot.Color = plot_color;
fault_slip_prob_plot.Units = 'centimeters';
fault_slip_prob_plot.Position = [6 8 11 9];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        'LineWidth',0.1,'Color',[0 0 0])
    hold on
end

% plotting seismicity data
M3=scatter(event_M3_coord_x/1000 , event_M3_coord_y/1000,7,'o','filled',...
    'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.5);
hold on
m45=scatter(event_M4_5_coord_x/1000 , event_M4_5_coord_y/1000,40,'pentagram','filled',...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],'LineWidth',0.5);


% Fault positions colorcoded with the slip probability at time t
if isempty(max_color_probability)==1
    color_idx = linspace(0,1,length(cmap_fault));                % index for assigning colors to each curve
else
    color_idx = linspace(0,max_color_probability,length(cmap_fault));
end

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

for i = 1:nr_fault

    [d , idx] = min(abs(color_idx -slip_prob_each_t(i,time_step)));

    plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',1.5,...
        'Color',cmap_fault(idx,:))

    hold on
end

% adding boxes of regions of interest for moderate magnitude earthquakes
% Fairview
plot([57 82 82 57 57], [154 154 174 174 154],'LineWidth',1,'Color',[0.64 0.078 0.18])
text(51.43,179.488,'d','Color',[0.64 0.078 0.18]);

% Pawnee
plot([218 243 243 218 218], [147 147 167 167 147],'LineWidth',1,'Color',[0.64 0.078 0.18])
text(234.035,174.16,'e','Color',[0.64 0.078 0.18]);

% Prague
plot([234 259 259 234 234], [47 47 67 67 47],'LineWidth',1,'Color',[0.64 0.078 0.18])
text(225.09,61.022,'b','Color',[0.64 0.078 0.18]);

% Cherokee
plot([85 110 110 85 85], [173 173 193 193 173],'LineWidth',1,'Color',[0.64 0.078 0.18])
text(85.74,199.94,'c','Color',[0.64 0.078 0.18]);

% Cushing
plot([230 255 255 230 230], [101 101 121 121 101],'LineWidth',1,'Color',[0.64 0.078 0.18])
text(222.105,105.13,'f','Color',[0.64 0.078 0.18]);

hold off

box('on');

cbar = colorbar('Position',[0.83 0.10 0.0381 0.83], ...
    'Ticks',[0 0.05 0.1 0.15 0.2 0.25 0.3],... 
    'LineWidth',0.2,...
    'FontName', 'Arial', ...
    'TickLabelInterpreter','latex',...
    'FontSize',colorbar_size);

set(get(cbar,'ylabel'),'string','Fault slip probability',...
    'Rotation',90,...
    'FontName', 'Arial', ...
    'interpreter', 'latex',...
    'fontsize',colorbar_size);

% set(gca,'ColorScale','log')


if isempty(max_color_probability)==1
    clim([0 1]);
else
    clim([0 max_color_probability]);
end

colormap(cmap_fault)

legend ([M3 m45], {['{\it M}' ' $\geq$ 3'] ['{\it M}' ' $\geq$ 4.5']},...
    'Location','northwest','Color',[1 1 1],'box','on'); 

axis('equal')

title('Dec. 2017',...
    'interpreter', 'latex')
ylabel('Y, Northting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);

axis([0 model_width 0 model_height])

text(216, 210,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,'Kansas', 'Interpreter','latex','FontSize', 8)
text(-39, 315,'a', 'Interpreter','latex','FontSize', 10)

axis([0 model_width 0 model_height])

xticks([0 100 200 300])
yticks([0 100 200 300])

set(gca,'Position',[0.04 0.10 0.8 0.83],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)

disp('spatial distribution of the faults colored by slip probability at the specified time is done')


%%%%%%%%%%%%%%%    Slip probability (Prague) - Fig. 6b     %%%%%%%%%%%%%%%%
Slip_probability_Prague = figure;
Slip_probability_Prague.Color = plot_color;
Slip_probability_Prague.Units = 'centimeters';
Slip_probability_Prague.Position = [6 8 5.4 4.2];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        'LineWidth',0.1,'Color',[0 0 0])
    hold on
end


% Fault positions colorcoded with the slip probability at time t
if isempty(max_color_probability)==1
    color_idx = linspace(0,1,length(cmap_fault));                % index for assigning colors to each curve
else
    color_idx = linspace(0,max_color_probability,length(cmap_fault));
end

M3=scatter(event_M3_coord_x/1000 , event_M3_coord_y/1000,7,'o','filled',...
    'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.5);
hold on

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

for i = 1:nr_fault
    [d , idx] = min(abs(color_idx -slip_prob_each_t(i,time_step)));
    
    plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',1.5,'Color',cmap_fault(idx,:))
    hold on
end

% plotting seismicity data
m45=scatter(event_M4_5_coord_x/1000 , event_M4_5_coord_y/1000,40,'pentagram','filled',...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],'LineWidth',0.5);
 
hold off

box('on');

if isempty(max_color_probability)==1
    clim([0 1]);
else
    clim([0 max_color_probability]);
end

colormap(cmap_fault)

axis('equal')

ylabel('Y, Northting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);

text(216, 210,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,'Kansas', 'Interpreter','latex','FontSize', 8)
text(-39, 315,'a', 'Interpreter','latex','FontSize', 10)


% % Aftershock
annotation(Slip_probability_Prague,'arrow',[0.515586419753087 0.46462962962963],...
    [0.64809722222222 0.640277777777778],'HeadWidth',6,'HeadLength',5);
% % Parague
annotation(Slip_probability_Prague,'arrow',[0.621419753086421 0.688055555555556],...
    [0.632978174603172 0.594920634920635],'HeadWidth',6,'HeadLength',5);
% % Foreshock
annotation(Slip_probability_Prague,'arrow',[0.629259259259259 0.699814814814815],...
    [0.736031746031746 0.756190476190475],'HeadWidth',6,'HeadLength',5);


% Prague 
text('String',['Prague',sprintf('\n'),'{\it M} = 5.8'],...
    'Position',[250.935828877006 56.7179144385027 0],'FontSize',7);
% Aftershock
text('String',['Aftershock'],...
    'Position',[234.55 58.376 0],'FontSize',7);
% Foreshock
text('String',['Foreshock'],...
    'Position',[251.508021390374 61.4508663101604 0],'FontSize',7);
% Wilzetta fault
text('String',['Wilzetta',sprintf('\n'),'fault'],...
    'Position',[246.905080213904 48.7967914438503 0],'FontSize',7 , 'Rotation',69);

% Parague header
text('String',['Parague 2011'],...
    'Position',[234.666 65.50 0],'FontSize',9);

text(229, 66.8,'b', 'Interpreter','latex','FontSize', 10)

axis([234 259 47 67])

xticks([240 250 260])
yticks([50 60])

set(gca,'Position',[0.20 0.127 0.73 0.92],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)


%%%%%%%%%%%%%      Slip probability (Cherokee) - Fig. 6c      %%%%%%%%%%%%%

Slip_probability_Cherokee = figure;
Slip_probability_Cherokee.Color = plot_color;
Slip_probability_Cherokee.Units = 'centimeters';
Slip_probability_Cherokee.Position = [6 8 5.4 4.2];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        'LineWidth',0.1,'Color',[0 0 0])
    hold on
end

% plotting seismicity data
M3=scatter(event_M3_coord_x/1000 , event_M3_coord_y/1000,7,'o','filled',...
    'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.5);
hold on
m45=scatter(event_M4_5_coord_x/1000 , event_M4_5_coord_y/1000,40,'pentagram','filled',...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],'LineWidth',0.5);


% Fault positions colorcoded with the slip probability at time t
if isempty(max_color_probability)==1
    color_idx = linspace(0,1,length(cmap_fault));                % index for assigning colors to each curve
else
    color_idx = linspace(0,max_color_probability,length(cmap_fault));
end

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

for i = 1:nr_fault
    [d , idx] = min(abs(color_idx -slip_prob_each_t(i,time_step)));
    
    plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',1.5,'Color',cmap_fault(idx,:))
    hold on
end
 
hold off

box('on');

if isempty(max_color_probability)==1
    clim([0 1]);
else
    clim([0 max_color_probability]);
end

colormap(cmap_fault)

axis('equal')

ylabel('Y, Northting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);

text(216, 210,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,'Kansas', 'Interpreter','latex','FontSize', 8)

% % Cherokee
annotation(Slip_probability_Cherokee,'arrow',[0.448950617283951 0.507746913580247],...
    [0.599960317460318 0.554603174603175],'HeadWidth',6,'HeadLength',5);

% Cherokee 
text('String',['Cherokee',sprintf('\n'),'{\it M} = 4.7'],...
    'Position',[96.0695187165775 181.598930481283 0],'FontSize',7);

% Cherokee header
text('String',['Cherokee 2015'],...
    'Position',[86 191.5 0],'FontSize',9);

text(79, 193,'c', 'Interpreter','latex','FontSize', 10)

axis([85 110 173 193])

xticks([90 100 110])
yticks([180 190])

set(gca,'Position',[0.20 0.127 0.73 0.92],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)


%%%%%%%%%%%%       Slip probability (Fairview) - Fig. 6d        %%%%%%%%%%%%

Slip_probability_Fairview = figure;
Slip_probability_Fairview.Color = plot_color;
Slip_probability_Fairview.Units = 'centimeters';
Slip_probability_Fairview.Position = [6 8 5.4 4.2];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        'LineWidth',0.1,'Color',[0 0 0])
    hold on
end


% plotting seismicity data
M3=scatter(event_M3_coord_x/1000 , event_M3_coord_y/1000,7,'o','filled',...
    'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.5);
hold on

% % unmapped fault
annotation(Slip_probability_Fairview,'line',[0.476388888888889 0.699814814814815],...
    [0.514285714285713 0.821706349206348],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'LineWidth',1.5,...
    'LineStyle','--');

hold on

m45=scatter(event_M4_5_coord_x/1000 , event_M4_5_coord_y/1000,40,'pentagram','filled',...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],'LineWidth',0.5);


% Fault positions colorcoded with the slip probability at time t
if isempty(max_color_probability)==1
    color_idx = linspace(0,1,length(cmap_fault));                % index for assigning colors to each curve
else
    color_idx = linspace(0,max_color_probability,length(cmap_fault));
end

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

% adding transparency to lines for better visualization
% transparency factor between 0 and 1; more transparent with smaller values
cmap_fault_trans = [cmap_fault(1:50,:) 0.85*ones(50,1);
              cmap_fault(51:256,:) ones(206,1)];

for i = 1:nr_fault
    [d , idx] = min(abs(color_idx -slip_prob_each_t(i,time_step)));
    
    plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',1.5,'Color',cmap_fault_trans(idx,:))
    hold on
end
 
hold off

box('on');

if isempty(max_color_probability)==1
    clim([0 1]);
else
    clim([0 max_color_probability]);
end

colormap(cmap_fault)

axis('equal')

ylabel('Y, Northting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);

text(216, 210,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,'Kansas', 'Interpreter','latex','FontSize', 8)
text(-39, 315,'a', 'Interpreter','latex','FontSize', 10)

% % Fairview foreshock
annotation(Slip_probability_Fairview,'arrow',[0.539020061728395 0.480308641975309],...
    [0.664960317460317 0.705793650793651],'HeadWidth',6,'HeadLength',5);

% % Fairview
annotation(Slip_probability_Fairview,'arrow',[0.624 0.679],...
    [0.596 0.535],'HeadWidth',6,'HeadLength',5);

% Fairview foreshcok
text('String',['Fairview',sprintf('\n'),'foreshock',sprintf('\n'),'{\it M} = 4.7'],...
    'Position',[59.7628877005347 167.196203208556 0],'FontSize',7);
% Fairview
text('String',['Fairview',sprintf('\n'),'{\it M} = 5.1'],...
    'Position',[73.93 161.66 0],'FontSize',7);

% Fairview header
text('String',['Fairview 2016'],...
    'Position',[57.5 172.5 0],'FontSize',9);

text(51, 173.5,'d', 'Interpreter','latex','FontSize', 10)

axis([57 82 154 174])

xticks([60 70 80])
yticks([160 170])

set(gca,'Position',[0.20 0.127 0.73 0.92],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)


%%%%%%%%%%%%%      Slip probability (Pawnee) - Fig. 6e        %%%%%%%%%%%%%

Slip_probability_Pawnee = figure;
Slip_probability_Pawnee.Color = plot_color;
Slip_probability_Pawnee.Units = 'centimeters';
Slip_probability_Pawnee.Position = [6 8 5.4 4.2];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        'LineWidth',0.1,'Color',[0 0 0])
    hold on
end

% plotting seismicity data
M3=scatter(event_M3_coord_x/1000 , event_M3_coord_y/1000,7,'o','filled',...
    'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.5);
hold on
m45=scatter(event_M4_5_coord_x/1000 , event_M4_5_coord_y/1000,40,'pentagram','filled',...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],'LineWidth',0.5);


% Fault positions colorcoded with the slip probability at time t
if isempty(max_color_probability)==1
    color_idx = linspace(0,1,length(cmap_fault));                % index for assigning colors to each curve
else
    color_idx = linspace(0,max_color_probability,length(cmap_fault));
end

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

for i = 1:nr_fault
    [d , idx] = min(abs(color_idx -slip_prob_each_t(i,time_step)));
    
    plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',1.5,'Color',cmap_fault(idx,:))
    hold on
end
 
hold off

box('on');

if isempty(max_color_probability)==1
    clim([0 1]);
else
    clim([0 max_color_probability]);
end

colormap(cmap_fault)

axis('equal')

ylabel('Y, Northting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);

text(216, 210,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,'Kansas', 'Interpreter','latex','FontSize', 8)
text(-39, 315,'a', 'Interpreter','latex','FontSize', 10)


% % Fairview
annotation(Slip_probability_Pawnee,'arrow',[0.558703703703704 0.61358024691358],...
    [0.579801587301588 0.458849206349207],'HeadWidth',6,'HeadLength',5);
% % unmapped fault
annotation(Slip_probability_Pawnee,'line',[0.5732421875 0.7841796875],...
    [0.60379797979798 0.529040404040404],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'LineWidth',1.5,...
    'LineStyle','--');

% Pawnee 
text('String',['Pawnee',sprintf('\n'),'{\it M} = 5.8'],...
    'Position',[231.303475935829 151.236631016042 0],'FontSize',7);
% Labette fault
text('String',['Labette fault'],...
    'Position',[237.352941176471 158.295454545454 0],'FontSize',7 , 'Rotation',62);
% Labette fault
text('String',['Stillwater fault'],...
    'Position',[222.205882352941 148.549465240642 0],'FontSize',7 , 'Rotation',57);

% Pawnee header
text('String',['Pawnee 2016'],...
    'Position',[218.8 165.3 0],'FontSize',9);


text(212.2, 167,'e', 'Interpreter','latex','FontSize', 10)

axis([218 243 147 167])

xticks([220 230 240])
yticks([150 160])

set(gca,'Position',[0.20 0.127 0.73 0.92],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)


%%%%%%%%%%%%%%      Slip probability (Cushing) - Fig. 6f      %%%%%%%%%%%%%

Slip_probability_Cushing = figure;
Slip_probability_Cushing.Color = plot_color;
Slip_probability_Cushing.Units = 'centimeters';
Slip_probability_Cushing.Position = [6 8 5.4 4.2];    % [x y w h]

axes('LineWidth',0.2);

time_series= dt:dt:time_project;
[~,time_step] = min(abs(time_series-time));

% county boarders (x data on the first half of columns
%                  y data on the second half of the columns)
for m=1:width(County_boarders)/2
    plot(County_boarders(:,m)/1000,County_boarders(:,m+width(County_boarders)/2)/1000,...
        'LineWidth',0.1,'Color',[0 0 0])
    hold on
end

% plotting seismicity data
M3=scatter(event_M3_coord_x/1000 , event_M3_coord_y/1000,7,'o','filled',...
    'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.5);
hold on
m45=scatter(event_M4_5_coord_x/1000 , event_M4_5_coord_y/1000,40,'pentagram','filled',...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],'LineWidth',0.5);


% Fault positions colorcoded with the slip probability at time t
if isempty(max_color_probability)==1
    color_idx = linspace(0,1,length(cmap_fault));                % index for assigning colors to each curve
else
    color_idx = linspace(0,max_color_probability,length(cmap_fault));
end

fault_points_x(:,1) = (fault_coord_x + fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_x(:,2) = (fault_coord_x - fault_length/2.*sin(deg2rad(fault_azi)))/1000;
fault_points_y(:,1) = (fault_coord_y + fault_length/2.*cos(deg2rad(fault_azi)))/1000;
fault_points_y(:,2) = (fault_coord_y - fault_length/2.*cos(deg2rad(fault_azi)))/1000;

for i = 1:nr_fault
    [d , idx] = min(abs(color_idx -slip_prob_each_t(i,time_step)));
    
    plot(fault_points_x(i,:) , fault_points_y(i,:),'LineWidth',1.5,'Color',cmap_fault(idx,:))
    hold on
end
 
hold off

box('on');

if isempty(max_color_probability)==1
    clim([0 1]);
else
    clim([0 max_color_probability]);
end

colormap(cmap_fault)

axis('equal')

ylabel('Y, Northting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);
xlabel('X, Easting [km]', ...
    'Interpreter', 'latex', ...
    'FontSize', axis_title_size);

text(216, 210,'Oklahoma', 'Interpreter','latex','FontSize', 8)
text(216, 229,'Kansas', 'Interpreter','latex','FontSize', 8)

% % Cushing
annotation(Slip_probability_Cushing,'arrow',[0.566543209876544 0.503827160493827],...
    [0.554603174603174 0.685634920634921],'HeadWidth',6,'HeadLength',5);
% % unmapped fault
annotation(Slip_probability_Cushing,'line',[0.652777777777778 0.503827160493827],...
    [0.56468253968254 0.463888888888889],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'LineWidth',1.5,...
    'LineStyle','--');

% Cushing 
text('String',['Cushing',sprintf('\n'),'{\it M} = 5'],...
    'Position',[237.560160427807 115.467914438503 0],'FontSize',7);

% Cushing header
text('String',['Cushing 2016'],...
    'Position',[231 119.4 0],'FontSize',9);

text(224.164438502674,120.6,'f', 'Interpreter','latex','FontSize', 10)

axis([230 255 101 121])

xticks([230 240 250])
yticks([110 120])

set(gca,'Position',[0.20 0.127 0.73 0.92],'FontSize',axis_size,...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1)


%% Saving figures
if save_figures=='Y'
    currentfolder = pwd;

    if isfolder([currentfolder,'\plots'])==0
        mkdir(currentfolder,'\plots');
    end

    filename=[currentfolder,'\plots\','deltap_contour'];
    print(deltap_contour, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(deltap_contour,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));
    
    filename=[currentfolder,'\plots\','fault_slip_prob_plot'];
    print(fault_slip_prob_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(fault_slip_prob_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));
    
    filename=[currentfolder,'\plots\','Histogram_plot'];
    print(Histogram_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(Histogram_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));

    filename=[currentfolder,'\plots\','cdf_plot'];
    print(cdf_plot, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(cdf_plot,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));

    filename=[currentfolder,'\plots\','cdf_plot_trunc'];
    print(cdf_plot_trunc, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(cdf_plot_trunc,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'));

    filename=[currentfolder,'\plots\','well_distribution'];
    print(well_distribution, '-dmeta', sprintf('-r%d', resolution), strcat(filename, '.emf'));
    print(well_distribution,'-dpdf', sprintf('-r%d', resolution), strcat(filename, '.pdf'))
end






