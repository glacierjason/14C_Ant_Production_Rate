%% Introduction


%% Load Data


% Clear workspace
clear, clc, close all

% Load supporting data
load supporting_data.mat

% Specify variables, sheet and range
opts = spreadsheetImportOptions("NumVariables", 11);
opts.Sheet = "Recalculated Ages";
opts.DataRange = "A3:K32";

% Specify column names and types
opts.VariableNames = ["Name", "elev", "compage", "compintunc", "compextunc", "tuage", "tuintunc", "tuextunc", "pubage", "pubintunc", "pubextunc"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double",  "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Name", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Name", "EmptyFieldRule", "auto");

% Find the data file
dataPath = fullfile(pwd, "PR_Data.xlsx");

% Load the data
Agesrecalcforlab = readtable(dataPath, opts, "UseExcel", false);

% Separate out variables for plotting
for i=1:width(Agesrecalcforlab)
    name = string(opts.VariableNames(i));
    ages.(name) = Agesrecalcforlab{:,i};
end

% Clear unnecessary data
clear opts i name Agesrecalcforlab dataPath

%% Figure 3

FieldNames = ["compage", "compextunc", "pubage", "pubextunc"];

color_map = containers.Map( ...
    {'pubage', 'tuage', 'compage',}, ...
    { [0, 0, 0], ... % Black
      [0.850, 0.372, 0.007], ... % Orange
      [0, 0, 0.8]});   % Dark Blue


% Figure settings
set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureWindowStyle','normal');

% Initialize
figure('Units','pixels','Position',[100 100 1000 1000],'Color','w');

% Set subplot locations manually for better plotting
fig1.left = 0.11; fig1.width = 0.75; fig1.h = 0.165;
fig1.b1 = 0.98 - 0.43;
fig1.b2 = fig1.b1 - (fig1.h + 0);
fig1.b3 = fig1.b2 - (fig1.h - 0.1);
fig1.b4 = fig1.b3 - (fig1.h - 0.02);
fig1.b5 = fig1.b4 - (fig1.h - 0.08);

ax = struct();
% Initialize axes based on subplot locations
ax.ax1 = axes('Position',[fig1.left fig1.b1 fig1.width 0.38], 'Box','off', 'Layer','top');
ax.ax2 = axes('Position',[fig1.left fig1.b2 fig1.width fig1.h], 'Box','off', 'Layer','top');
ax.ax3 = axes('Position',[fig1.left fig1.b3 fig1.width fig1.h], 'Box','off', 'Layer','top');
ax.ax4 = axes('Position',[fig1.left fig1.b4 fig1.width fig1.h], 'Box','off', 'Layer','top');
ax.ax5 = axes('Position',[fig1.left fig1.b5 fig1.width fig1.h], 'Box','off', 'Layer','top');

% Create and array of axes and hold on
axes_arr = [ax.ax1 ax.ax2 ax.ax3 ax.ax4 ax.ax5];
hold(axes_arr, 'on')


% - ax1 -
ages.tam_subset_samples = [1, 2, 7, 8, 9, 10];
agecolumns = [1,3];

for j = agecolumns
    samples = ages.(FieldNames(j));
    uncs = ages.(FieldNames(j+1));
    color = color_map(FieldNames(j));

    for i = ages.tam_subset_samples
        errorbar(ax.ax1, samples(i), ages.elev(i), ...
            uncs(i), uncs(i), ...
            'horizontal','LineStyle','none', 'LineWidth', 1, 'Color', color, 'CapSize',10);
        scatter(ax.ax1, samples(i),  ages.elev(i), 160, color, 'filled', 'Marker','o', 'MarkerFaceAlpha',1);
    end
end


% Plot the pchip spline for the TAMs
splines.tam_spline.x_comp = [sort(ages.compage(ages.tam_subset_samples),'descend')];
splines.tam_spline.y = [sort(ages.elev(ages.tam_subset_samples), 'descend')];
splines.tam_spline.pp_comp = pchip(splines.tam_spline.x_comp, splines.tam_spline.y);
splines.tam_spline.xx_comp = linspace(4200, 15700, 10000);
plot(ax.ax1, splines.tam_spline.xx_comp, ppval(splines.tam_spline.pp_comp, splines.tam_spline.xx_comp), 'LineWidth', 2, 'Color','b');
splines.tam_spline.yy_comp = ppval(splines.tam_spline.pp_comp, splines.tam_spline.xx_comp);


for i = ages.tam_subset_samples
    errorbar(ax.ax1, ages.pubage(i), ages.elev(i), ...
        ages.pubextunc(i), ages.pubextunc(i), ...
        'horizontal','LineStyle','none', 'LineWidth', 1, 'Color', 'k', 'CapSize',10);
    scatter(ax.ax1, ages.pubage(i),  ages.elev(i), 160, 'k', 'filled', 'Marker','o', 'MarkerFaceAlpha',1);
end

splines.tam_spline.x_orig = [sort(ages.pubage(ages.tam_subset_samples),'descend')];
splines.tam_spline.pp_orig = pchip(splines.tam_spline.x_orig, splines.tam_spline.y);
splines.tam_spline.xx_orig = linspace(5100, 60000, 10000);
plot(ax.ax1, splines.tam_spline.xx_orig, ppval(splines.tam_spline.pp_orig, splines.tam_spline.xx_orig), 'LineWidth', 2, 'Color', [0 0 0 0.8]);
splines.tam_spline.yy_orig = ppval(splines.tam_spline.pp_orig, splines.tam_spline.xx_orig);

uistack(fill(ax.ax1, [splines.tam_spline.xx_comp(:); flipud(splines.tam_spline.xx_orig(:))], ...
    [splines.tam_spline.yy_comp(:); flipud(splines.tam_spline.yy_orig(:))], ...
    'b', 'FaceAlpha', 0.05, 'EdgeColor', 'none'), 'top');

scatter(ax.ax1, 24650, 558, 160, 'k', 'filled', 'Marker','>')


t_rate = struct();
t_rate.a = linspace(18000, 21000, 2);
t_rate.b = linspace(18000, 20500, 2);
t_rate.c = linspace(18000, 19000, 2);
t_rate.d = linspace(18000, 18200, 2); 
rates_cm_per_yr = [1, 2, 10, 50];
rates_m_per_yr  = rates_cm_per_yr / 100;
y0 = 100;  

plot(ax.ax1, t_rate.a, (y0 + (t_rate.a - t_rate.a(1)) * rates_m_per_yr(1)), ...
    'Color', 'k', 'LineWidth', 1.8);
plot(ax.ax1, t_rate.b, (y0 + (t_rate.b - t_rate.b(1)) * rates_m_per_yr(2)), ...
    'Color', 'k', 'LineWidth', 1.8);
plot(ax.ax1, t_rate.c, (y0 + (t_rate.c - t_rate.c(1)) * rates_m_per_yr(3)), ...
    'Color', 'k', 'LineWidth', 1.8);
plot(ax.ax1, t_rate.d, (y0 + (t_rate.d - t_rate.d(1)) * rates_m_per_yr(4)), ...
    'Color', 'k', 'LineWidth', 1.8);

temp.unique_categories = {'Compilation'};
temp.legend_colors = [ ...
     0, 0, 0.8]; % Compilation

% Initialize dummy patch handles for legend
temp.legend_handles = gobjects(1, length(temp.unique_categories));

% Create patches for legend
for i = 1:length(temp.unique_categories)
    temp.legend_handles(i) = plot(ax.ax1, NaN, NaN, 'o', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', temp.legend_colors(i,:), ...
        'MarkerEdgeColor', 'none');
end

temp.unique_categories2 = {'Published'};

% Initialize dummy patch handles for published data
temp.legend_handles2 = gobjects(1, 1);

% Create patches for legend
temp.legend_handles2 = plot(ax.ax1, NaN, NaN, 'o', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'none');

% Generate the legend with transparent patches
leg = legend(ax.ax1, [temp.legend_handles, temp.legend_handles2], ...
    [temp.unique_categories, temp.unique_categories2], ...
   'Location', [0.708095100563742,0.772556867923571,0.146789665160481,0.074820491373734]);

leg.ItemTokenSize = [15, 15]; % Set size for legend markers (width, height)

text(ax.ax1,15500,230,'50 cm/yr','FontSize',18,'FontName','Helvetica','Clipping','off');
text(ax.ax1, 19000, 225, '10 cm/yr','FontSize',18,'FontName','Helvetica','Clipping','off');
text(ax.ax1,20600,175,'2 cm/yr','FontSize',18,'FontName','Helvetica','Clipping','off');
text(ax.ax1, 21100, 135, '1 cm/yr','FontSize',18,'FontName','Helvetica','Clipping','off');

clear temp

ylabel(ax.ax1,'Elevation (m asl)');
ylim(ax.ax1, [0, 600])


% - ax2 -
cla(ax.ax2);
plot(ax.ax2, clim_data.waisO18.age, clim_data.waisO18.oxy, 'k-', 'LineWidth',3, 'DisplayName', 'WAIS Divide \delta^{18}O (Steig et al., 2015)');
ylabel(ax.ax2,'\delta^{18}O (‰)','Interpreter','tex');
set(ax.ax2,'YAxisLocation','left','YLim',[-42.3 -32], ...
        'TickDir','out','LineWidth',1.25,'FontName','Helvetica', ...
        'Color','none','PositionConstraint','innerposition');
legend(ax.ax2, 'location', [0.443863749254006,0.519678518103771,0.408752237982982,0.037096240241899])


% - ax3 -
ax.ax3_1 = plot(ax.ax3, clim_data.EATA.age, clim_data.EATA.TA, 'Color',[0 0.7 0.5], ...
    'LineWidth',3, 'DisplayName', 'East Antarctica Temp. (Parrenin et al., 2013)');
ylim(ax.ax3, [-10 2]);
ax.ax3.YColor = [0 0.7 0.5];
ax.ax3.YLabel.Color = [0 0.7 0.5];
ylabel(ax.ax3,'Temp Anomaly (°C)'); yticks([-10 -8 -6 -4 -2 0 1]);
legend(ax.ax3, [ax.ax3_1], 'location', [0.118222497681137,0.352071977373951,0.476651100733822,0.032089552238806])


% - ax4 -
ax.ax4_1 = plot(ax.ax4, clim_data.creelSLR.age, clim_data.creelSLR.mid, ...
    'Color',[0 0.2 0.7], 'LineWidth', 3, 'DisplayName', 'Creel et al., 2024');
fill(ax.ax4, [clim_data.creelSLR.age; flipud(clim_data.creelSLR.age)], ...
    [clim_data.creelSLR.lowunc; flipud(clim_data.creelSLR.highunc)], [0 0.2 0.7], 'FaceAlpha', 0.3)
ax.ax4_2 = plot(ax.ax4, clim_data.lambeckSLR.age, clim_data.lambeckSLR.SLR, ...
    'Color',[0 0.2 0.7], 'LineWidth', 3, 'LineStyle', '--','DisplayName', 'Lambeck et al., 2014');
ylabel(ax.ax4,'GMSL (m)');
ylim(ax.ax4,[-150 2]);
ax.ax4.YColor = [0 0.2 0.7];
ax.ax4.YLabel.Color = [0 0.2 0.7];
yticks([-150 -100 -50 0])
legend(ax.ax4, [ax.ax4_1 ax.ax4_2], 'location', [0.601271911003733,0.226834282807201,0.251834556370914,0.060261194029851])


% - ax5 -
ax.ax5_1 = plot(ax.ax5, clim_data.solin65n.age, clim_data.solin65n.insolref, 'LineWidth',3, ...
    'Color', [0.9 0.5 0.1], 'LineStyle', '--', 'DisplayName', '65˚N Summer');
ax.ax5_2 = plot(ax.ax5, clim_data.solin65s.age, clim_data.solin65s.insolref, 'LineWidth',3, ...
    'Color', [0.9 0.5 0.1], 'DisplayName', '65˚S Summer');
ax.ax5_3 = plot(ax.ax5, clim_data.solinavg.age, clim_data.solinavg.insolref, 'k--', 'LineWidth', 2, ...
    'DisplayName', 'Whole Earth Avg.');
ylabel(ax.ax5,'Insolation (% rel. to present)', 'HorizontalAlignment', 'center');
ylim(ax.ax5,[0.95 1.15]);
yticks([0.95 1 1.05 1.1 1.15])
ax.ax5.YColor = [0.9 0.5 0.1];
ax.ax5.YLabel.Color = [0.9 0.5 0.1];
xticks([0 5000 10000 15000 20000 25000])
ax.ax5.XAxis.Exponent = 0; xlabel(ax.ax5,'Time (yr BP)');
legend(ax.ax5, [ax.ax5_1, ax.ax5_2, ax.ax5_3], 'location', [0.117306223835577,0.21095275366838,0.219146097398265,0.088432835820895])


% Linkaxes and x-axis formatting
linkaxes(axes_arr,'x');
set([ax.ax1 ax.ax2 ax.ax3 ax.ax4 ax.ax5], 'FontSize',18, 'FontName','Helvetica', ...
    'TickDir','out', 'LineWidth',1.25, 'Color','none');
set([ax.ax1 ax.ax2 ax.ax3 ax.ax4], 'XColor','none','XTick',[]);
xlim([0 25000]);

% Set y-axis sides to alternate for spacing
set([ax.ax1 ax.ax3 ax.ax5],'YAxisLocation','right');
set([ax.ax2 ax.ax4],'YAxisLocation','left');

% y-label spacing to solve axis labeling overlap issues
for ax100 = [ax.ax2 ax.ax4]
    ax100.YLabel.Units = 'normalized'; fig1.p = ax100.YLabel.Position; fig1.p(1) = -0.087; 
    ax100.YLabel.Position = fig1.p;
end
for ax100 = [ax.ax1 ax.ax3]
    ax100.YLabel.Units = 'normalized'; fig1.p = ax100.YLabel.Position; fig1.p(1) = 1.07;  ax100.YLabel.Position = fig1.p;
end
for ax100 = [ax.ax5]
    ax100.YLabel.Units = 'normalized'; fig1.p = ax100.YLabel.Position; fig1.p(1) = 1.11; 
    fig1.p(2) = 0.6; ax100.YLabel.Position = fig1.p;
end


% Set axis for multipanel overlay data
fig1.Ytop = ax.ax1.Position(2) + ax.ax1.Position(4);   % top of top panel
fig1.Ybot = ax.ax5.Position(2);                     % bottom of lowest panel
bandAx = axes('Position',[fig1.left fig1.Ybot fig1.width (fig1.Ytop - fig1.Ybot)], ...
    'Color','none','HitTest','off','HandleVisibility','off', ...
    'XLim',[0 25000],'YLim',[0 1], 'XTick',[], 'YTick',[], ...
    'XColor','none','YColor','none','Layer','bottom');

% Assign the size of multipanel bands
% y0,y1 are fractions of the bandAx height (0 = bottom, 1 = top)
fig1.span_all = [0 1];
fig1.bands = [4000  7000;
         13700 14200];
fig1.cols  = [1 0 0; 1 0 0];
fig1.alpha_main = 0.30;

for b = 1:size(fig1.bands,1)
    fig1.x0 = fig1.bands(b,1); fig1.x1 = fig1.bands(b,2);
    fig1.y0y1 = fig1.span_all;
    patch(bandAx,[fig1.x0 fig1.x1 fig1.x1 fig1.x0],[fig1.y0y1(1) fig1.y0y1(1) fig1.y0y1(2) fig1.y0y1(2)], ...
          fig1.cols(b,:), 'FaceAlpha', fig1.alpha_main, 'EdgeColor','none');
end

% keep behind data and lock x with the other axes
linkaxes([axes_arr bandAx],'x');
% move the entire background axis behind all other axes
uistack(bandAx, 'bottom');

% labels
fig1.yl1 = ylim(ax.ax1);
text(ax.ax1,12600,fig1.yl1(2)*1.05,'MWP-1A','FontSize',18,'FontName','Helvetica','Clipping','off');
text(ax.ax1, 300,fig1.yl1(2)*1.05,'Holocene Thermal Maximum','FontSize',18,'FontName','Helvetica','Clipping','off');

fig = gcf;
% Save a high-resolution figure
%exportgraphics(fig, 'Figure3.jpg', 'Resolution', 300)
hold off;
clear fig fig1 ax ax100 axes_arr b i bandAx



%% 
figure('Units','pixels','Position',[100 100 1000 500],'Color','w');
plot_subset = [22:30];
landform_age = struct();
landform_age.blue_orig = [13400, 300]; % Ketelers
landform_age.blue_comp = [9600, 200];
landform_age.dry_orig = [22000, 1300];
landform_age.dry_comp = [13000, 400];
landform_age.utse_orig = [4400, 200];
landform_age.utse_comp = [3600, 200];

errorbar(ages.pubage(plot_subset), ages.elev(plot_subset), ...
    ages.pubintunc(plot_subset), ages.pubintunc(plot_subset), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.698, 0.094, 0.169], 'CapSize', 10, 'HandleVisibility', 'off');
hold on;
errorbar(ages.compage(plot_subset), ages.elev(plot_subset), ...
    ages.compintunc(plot_subset), ages.compintunc(plot_subset), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.129, 0.400, 0.675], 'CapSize', 10, 'HandleVisibility', 'off');
h1 = scatter(ages.pubage(plot_subset), ages.elev(plot_subset), 70, [0.698, 0.094, 0.169], 'filled', 'DisplayName', 'Published age');
h2 = scatter(ages.compage(plot_subset), ages.elev(plot_subset), 70, [0.129, 0.400, 0.675], 'filled', 'DisplayName', 'Updated age');

errorbar(landform_age.blue_orig(1), 1425, ...
    landform_age.blue_orig(2), landform_age.blue_orig(2), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.937, 0.541, 0.384], 'CapSize', 10, 'HandleVisibility', 'off');
h3 = scatter(landform_age.blue_orig(1), 1425, 100, [0.937, 0.541, 0.384], 'filled', 'v', 'DisplayName', 'Landform age (published)');
errorbar(landform_age.dry_orig(1), 1620, ...
    landform_age.dry_orig(2), landform_age.dry_orig(2), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.937, 0.541, 0.384], 'CapSize', 10, 'HandleVisibility', 'off');
scatter(landform_age.dry_orig(1), 1620, 100, [0.937, 0.541, 0.384], 'filled', 'v', 'HandleVisibility', 'off')
errorbar(landform_age.utse_orig(1), 1360, ...
    landform_age.utse_orig(2), landform_age.utse_orig(2), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.937, 0.541, 0.384], 'CapSize', 10, 'HandleVisibility', 'off');
scatter(landform_age.utse_orig(1), 1360, 100, [0.937, 0.541, 0.384], 'filled', 'v', 'HandleVisibility', 'off')

errorbar(landform_age.blue_comp(1), 1425, ...
    landform_age.blue_comp(2), landform_age.blue_comp(2), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.404, 0.663, 0.812], 'CapSize', 10, 'HandleVisibility', 'off');
h4 = scatter(landform_age.blue_comp(1), 1425, 100, [0.404, 0.663, 0.812], 'filled', 'v', 'DisplayName', 'Landform age (updated)');
errorbar(landform_age.dry_comp(1), 1620, ...
    landform_age.dry_comp(2), landform_age.dry_comp(2), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.404, 0.663, 0.812], 'CapSize', 10, 'HandleVisibility', 'off');
scatter(landform_age.dry_comp(1), 1620, 100, [0.404, 0.663, 0.812], 'filled', 'v', 'HandleVisibility', 'off')
errorbar(landform_age.utse_comp(1), 1360, ...
    landform_age.utse_comp(2), landform_age.utse_comp(2), ...
'horizontal','LineStyle','none', 'LineWidth', 2, 'Color', [0.404, 0.663, 0.812], 'CapSize', 10, 'HandleVisibility', 'off');
scatter(landform_age.utse_comp(1), 1360, 100, [0.404, 0.663, 0.812], 'filled', 'v', 'HandleVisibility', 'off')

scatter([24750, 24750, 24750, 24750], [1607, 1395, 1401, 1394], 70, [0.698, 0.094, 0.169], 'filled', 'Marker', '>')
scatter([24400, 24400], [1395, 1401], 70, [0.129, 0.400, 0.675], 'filled', 'Marker', '>')


text(1300, 1390, 'Utsteinen (\Delta 0.8 ka)', 'FontSize', 18, 'FontWeight', 'normal', 'FontName', 'Helvetica', 'Color', [0.4 0.4 0.4])
text(13300, 1455, 'Ketelers (\Delta 3.8 ka)', 'FontSize', 18, 'FontWeight', 'normal', 'FontName', 'Helvetica', 'Color', [0.4 0.4 0.4])
text(15300, 1650, 'Dry (\Delta 9.0 ka)', 'FontSize', 18, 'FontWeight', 'normal', 'FontName', 'Helvetica', 'Color', [0.4 0.4 0.4])


xlim([0, 25000])
ylim([1300, 1700])
xlabel('Time (yr BP)', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Helvetica');
ylabel('Elevation (m)', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Helvetica');
set(gca, 'FontSize', 18, 'FontName', 'Helvetica');
ax = gca;
xticks([0 5000 10000 15000 20000 25000])
ax.XAxis.Exponent = 0;
box on;

legend([h1, h2, h3, h4], 'Location', 'northwest', 'FontSize', 18, 'Box', 'on');

fig = gcf;
%exportgraphics(fig, 'Figure4.jpg', 'Resolution', 300)
hold off;



%% Thinning Rates

splines.tam_spline.pp_comp_deriv = fnder(splines.tam_spline.pp_comp);
splines.tam_spline.slope_comp = ppval(splines.tam_spline.pp_comp_deriv, ...
    splines.tam_spline.xx_comp) * 1000;

splines.tam_spline.pp_orig_deriv = fnder(splines.tam_spline.pp_orig);
splines.tam_spline.slope_orig = ppval(splines.tam_spline.pp_orig_deriv, ...
    splines.tam_spline.xx_orig) * 1000;


% Plotting
% Set figure settings
set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureWindowStyle','normal');   % not docked

figure('Color','w');


% TAM

hold on; box on;
a = plot(splines.tam_spline.xx_comp/1000, splines.tam_spline.slope_comp/10, ...
    'LineWidth',2, 'Color','b');
b = plot(splines.tam_spline.xx_orig/1000,  splines.tam_spline.slope_orig/10, ...
    '--', 'LineWidth',2, 'Color',[0 0 1 0.6]);

xlim([0 25]); ylim([0 20]);
ax = gca;
ax.XAxis.Exponent = 0;
xlabel('Age (kyr)');
ylabel('Thinning rate (cm/yr)')
%title('Thinning Rate');
legend([a b], 'Compilation PR Thinning Rate', 'Published Data Thinning Rate');
set(gca, 'FontSize', 24, 'FontName', 'Helvetica', 'Box', 'on','YGrid', 'off');
grid on;
set(gcf,'Position',[100 100 1200 800]);




%% Calculate and plot residuals

% Calculate residuals
residuals = ages.pubage(:) - ages.compage(:);                    
residual_unc = sqrt((ages.compextunc(:)).^2 + (ages.pubextunc(:)).^2 ); % propagate errors
residuals = residuals/1e4;
residual_unc = residual_unc/1e4;

% Create bar plot with error bars
figure('Color','w','Units','pixels','Position',[200 200 1200 500]);
hb = bar(residuals);
hold on;

% Add errorbars
x = hb.XEndPoints;
errorbar(x, residuals, residual_unc, 'k.', 'LineWidth', 1.2, 'CapSize', 10);

% Aesthetics
xlabel('Sample');
ylabel('Residual (published - calibrated) [kyr]');
grid on;
set(gca,'FontSize',14,'TickDir','out');

% Replace x-ticks with sample names (rotate for readability)
xticks(1:numel(residuals));
xticklabels(ages.Name);
xtickangle(45);

ax = gca;
ax.Box = 'off';               % remove the box (top/right border)
ax.XAxisLocation = 'bottom';  % ensure x ticks only on bottom
ax.YAxisLocation = 'left';    % ensure y ticks only on left

hold off;




%% Calculate and plot residuals

% Calculate residuals
residuals = ages.pubage(:) - ages.tuage(:);                    
residual_unc = sqrt((ages.tuextunc(:)).^2 + (ages.pubextunc(:)).^2 ); % propagate errors
residuals = residuals/1e4;
residual_unc = residual_unc/1e4;

% Create bar plot with error bars
figure('Color','w','Units','pixels','Position',[200 200 1200 500]);
hb = bar(residuals);
hold on;

% Add errorbars
x = hb.XEndPoints;
errorbar(x, residuals, residual_unc, 'k.', 'LineWidth', 1.2, 'CapSize', 10);

% Aesthetics
xlabel('Sample');
ylabel('Residual (published - calibrated) [kyr]');
grid on;
set(gca,'FontSize',14,'TickDir','out');

% Replace x-ticks with sample names (rotate for readability)
xticks(1:numel(residuals));
xticklabels(ages.Name);
xtickangle(45);

ax = gca;
ax.Box = 'off';               % remove the box (top/right border)
ax.XAxisLocation = 'bottom';  % ensure x ticks only on bottom
ax.YAxisLocation = 'left';    % ensure y ticks only on left

hold off;