% This code was originally written by Keir Nichols as a PhD student at
% Tulane University to calculate the spallation production rate of in-situ 14C in
% quartz using the CRONUS-A measurement as a calibration measurement. It was
% subsequently updated by Jason Drebber as a PhD Student at the Colorado
% School of Mines to include a larger number of CRONUS-A measurements from
% different labs and expanded to incorporate the LSDn scaling framework and
% more statistical tests

%% Set-Up
% turn this on to clear figures, output in command window and workspace variables
clf, clc, clear, close all

% Get the time it takes to run the whole code
tic;

% Use a seed to produce the same output from random methods everytime
rng(84735)


%% Define Constants and Load Data
const.constmax.lambda = log(2)/5670; % 14C decay constants (calculated using a half-life of
const.lambda = log(2)/5700;    % of 5700 +/- 30 years following Hippe and Lifton, 2014;
const.constmin.lambda = log(2)/5730; % This is the IAEA accepted value as of June 2024)

const.L = 267.8; % Effective e-folding length in atmospheric pressure (hPa) calculated for in-situ 14C (needed for muon production)
const.constmin.L = 267.8; % Used for the uncertainty analysis
const.constmax.L = 267.8; % Used for the uncertainty analysis

const.Fsp = 1; % Used in stone scaling to account for muogenic production (Balco, 2008)
const.constmax.Fsp = 1; % Used for the uncertainty analysis
const.constmin.Fsp = 1; % Used for the uncertainty analysis

CA.lat = -77.883; % Latitude of CRONUS-A sample (Jull et al., 2015)
CA.long = 160.9431; % Longitude of CRONUS-A sample (Jull et al., 2015)
CA.z = 1612; % Elevation of CRONUS-A site in meters (Jull et al., 2015)
%CA.z = 1679; % Elevation of CRONUS-A site in meters (Balco et al., 2019)
CA.z_error = 10; % Error in elevation measurements from a GPS collected point (Lecavlier, 2022)

const.atmp = ERA40atm(CA.lat, CA.long, CA.z); % Calculate the atmospheric pressure of the CRONUS-A site
const.constmin.atmp = ERA40atm(CA.lat, CA.long, CA.z+CA.z_error); % Used for the uncertainty analysis
const.constmax.atmp = ERA40atm(CA.lat, CA.long, CA.z-CA.z_error); % Used for the uncertainty analysis

% Load data
load CRONUSA.mat % CRONUS-A measurements
load extlab.mat % Extraction Lab data
CA.conc = CRONUSA(:,4); % CA is a structure of CRONUS-A data; conc in the concentration measurements
CA.error = CRONUSA(:,5); % Published error on the concentration measurements

% Load the data by extraction date
date_extracted.date = load("extractiondate.mat", 'DateExtracted');
date_extracted.conc = load("dateext_conc.mat");
date_extracted.lab = load('dateext_lab.mat');
date_extracted.AMS = load('dateext_AMS.mat');
date_extracted.date = date_extracted.date.DateExtracted;
date_extracted.conc = date_extracted.conc.dateext_conc;
date_extracted.lab = date_extracted.lab.dateext_lab;
date_extracted.AMS = date_extracted.AMS.dateext_AMS;


% Assign a color map (dictionary) that corresponds to each lab for plotting
color_map = containers.Map( ...
    {'Purdue', 'P-CEGS', 'TU-CEGS', 'LDEO', 'ANSTO-UOW', ...
     'ETH Zurich', 'Intercomparison', 'University of Cologne'}, ...
    { [0.121, 0.467, 0.706], ... % Blue
      [0.172, 0.627, 0.745], ... % Teal
      [0.850, 0.372, 0.007], ... % Orange
      [0.702, 0.871, 0.412], ... % Light Green
      [0.525, 0.396, 0.750], ... % Purple
      [0.984, 0.502, 0.447], ... % Coral
      [0.596, 0.596, 0.596], ... % Grey
      [0.255, 0.255, 0.255]});   % Dark Grey

disp('Data Loaded...')


%% Create a matrix of the lab data
% This is needed for plotting the data as a violin plot

%Initialize the matrix data
CA.matrix_data = createMatrix(extlab, CA.conc);

%Reorganize the matrix corresponding to the order of the plotting
temp.c1 = CA.matrix_data(:,6);
temp.c2 = CA.matrix_data(:,5);
temp.c3 = CA.matrix_data(:,7);
temp.c4 = CA.matrix_data(:,4);
temp.c5 = CA.matrix_data(:,1);
temp.c6 = CA.matrix_data(:,2);
temp.c7 = CA.matrix_data(:,3);
temp.c8 = CA.matrix_data(:,8);

CA.matrix_data = [temp.c1, temp.c2, temp.c3, temp.c4, ...
                  temp.c5, temp.c6, temp.c7, temp.c8];

clear temp;


%% Put the data into a structured array
% This will store the data for each lab in it's own array and make managing
% the data easier in the rest of the code

% Store the unique data for each lab in the lab structure
[lab, const] = structureLabData(extlab, CA.conc, CA.error, const);

% Store the lab names 
const.labs = {"Purdue", "P-CEGS", 'TU-CEGS', 'LDEO', 'ANST0-UOW', 'ETH-Zurich', 'Intercomparison', 'University of Cologne'};
const.labs_mod = {"Purdue", "P-CEGS", 'TU-CEGS', 'LDEO', 'ANST0-UOW', 'ETH-Zurich', 'Intercomp.', 'U. of Cologne'};

% Store lab names again for uncertainty analysis
const.constmin.labs = {"Purdue", "P-CEGS", 'TU-CEGS', 'LDEO', 'ANST0-UOW', 'ETH-Zurich', 'Intercomparison', 'University of Cologne'};
const.constmin.labs_mod = {"Purdue", "P-CEGS", 'TU-CEGS', 'LDEO', 'ANST0-UOW', 'ETH-Zurich', 'Intercomp.', 'U. of Cologne'};

% Store lab names again for uncertainty analysis
const.constmax.labs = {"Purdue", "P-CEGS", 'TU-CEGS', 'LDEO', 'ANST0-UOW', 'ETH-Zurich', 'Intercomparison', 'University of Cologne'};
const.constmax.labs_mod = {"Purdue", "P-CEGS", 'TU-CEGS', 'LDEO', 'ANST0-UOW', 'ETH-Zurich', 'Intercomp.', 'U. of Cologne'};


% Reorder the data to correspond to the lab name order
lab = lab([6, 5, 7, 4, 1, 2, 3, 8]);

% Store the lab data without the University of Cologne
labtrim = lab(1:7);

% Save the raw data including the University of Cologne for statistics
CA.conc_raw = CA.conc;
CA.error_raw = CA.error;
CA.extlab_raw = extlab;

% Remove the University of Cologne data from the compiled data because
% there is low certainty in the values
[valid_idx] = find(extlab ~= "University of Cologne");
CA.conc = CA.conc(valid_idx);
CA.error = CA.error(valid_idx);
extlab = extlab(valid_idx);

clear valid_idx


%% Plot the main compilation figure
% Figure 1
% This figure shows the compiled CRONUS-A kernel density for each lab in a
% single figure with each lab shown as a different color. The labs are
% stacked on top of each other for easy comparison.

% Define the range for plotting
temp.buffer = 1 * (max(CA.conc) - min(CA.conc)); % Add buffer on both sides so the figure is not cut off
conc_range.full_comp = linspace(min(CA.conc) - temp.buffer, max(CA.conc) + temp.buffer, 1000);

% Scale x-axis for display
conc_range.full_comp_scaled = conc_range.full_comp/1e5; % Scale range for display purposes

% Calculate kernel densities for each lab
for i = 1:length(labtrim)
    lab_name = labtrim(i).name;

    % Store kernel densities and metadata in the new variable
    labtrim(i).color = color_map(lab_name);
    labtrim(i).y = normpdf(conc_range.full_comp, labtrim(i).conc, labtrim(i).error);
    labtrim(i).totalconc = sum(labtrim(i).y, 1); % Sum across all densities
end

figure;
hold on;

% Initialize cumulative density
temp.tot_density = zeros(1, length(conc_range.full_comp));

% Initialize cumulative density
temp.tot_scale_density = zeros(1, length(conc_range.full_comp_scaled));

% Loop through each lab to compute and stack densities
for i = 1:length(labtrim)
    if i==3 % Pick the Tulane data and plot it differently for emphasis if needed
            % Add the current density to the cumulative density
            temp.stacked_density = temp.tot_density + labtrim(i).totalconc;

            % Add the current scaled density to the scaled cumulative density
            temp.stacked_scale_density = temp.tot_scale_density + labtrim(i).totalconc/1e5;

            % Plot the stacked area
            fill([conc_range.full_comp_scaled, fliplr(conc_range.full_comp_scaled)], ...
            [temp.stacked_scale_density, fliplr(temp.tot_scale_density)], ...
            labtrim(i).color, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

            % Update cumulative density
            temp.tot_density = temp.stacked_density;


            % Update scaled cumulative density
            temp.tot_scale_density = temp.stacked_scale_density;

    else
            % Add the current density to the cumulative density
            temp.stacked_density = temp.tot_density + labtrim(i).totalconc;

            % Add the current scaled density to the scaled cumulative density
            temp.stacked_scale_density = temp.tot_scale_density + labtrim(i).totalconc/1e5;

            % Plot the stacked area
            fill([conc_range.full_comp_scaled, fliplr(conc_range.full_comp_scaled)], ...
            [temp.stacked_scale_density, fliplr(temp.tot_scale_density)], ...
            labtrim(i).color, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

            % Update cumulative density
            temp.tot_density = temp.stacked_density;


            % Update scaled cumulative density
            temp.tot_scale_density = temp.stacked_scale_density;
    end
end

% Store the total density in the CA structure
CA.tot_density = temp.tot_density;

% Add axis labels with refined font sizes
xlabel('^{14}C concentration (10^{5} atoms g^{-1})', 'FontSize', 36, 'FontWeight', 'Bold', 'FontName', 'Helvetica');
ylabel('Probability Density', 'FontSize', 36, 'FontWeight', 'Bold', 'FontName', 'Helvetica');

ax = gca;  % Get current axes
ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
ax.XAxis.Exponent = 0;
ax.XAxis.TickLabelFormat = '%.1f'; % Keep scientific notation
ax.YTickLabel = [];
ax.YTick = [];

% Define legend categories and colors
temp.unique_categories = {'Intercomparison', 'ETH Zurich', 'ANSTO-UOW', ...
                          'LDEO', 'TU-CEGS', 'P-CEGS', 'Purdue'};

legend_colors = [ ...
     0.596, 0.596, 0.596;  % Intercomparison
     0.984, 0.502, 0.447;  % ETH Zurich
     0.525, 0.396, 0.750;  % ANSTO-UOW
     0.702, 0.871, 0.412;  % LDEO
     0.850, 0.372, 0.007;  % TU-CEGS
     0.172, 0.627, 0.745;  % P-CEGS
     0.121, 0.467, 0.706]; % Purdue


% Initialize dummy patch handles for legend
temp.legend_handles = gobjects(1, length(temp.unique_categories));

% Create patches for legend
for i = 1:length(temp.unique_categories)
    temp.legend_handles(i) = patch(NaN, NaN, legend_colors(i, :), ...
        'FaceAlpha', 0.6, 'EdgeColor', 'none');
end

% Generate the legend with transparent patches
leg = legend(temp.legend_handles, temp.unique_categories, ...
    'Location', 'northwest', ...
    'FontSize', 36, 'Box', 'off', 'TextColor', 'k', 'FontName', 'Helvetica');

% [left, bottom, width, height] in normalized units
leg.Position = [0.25, 0.50, 0.1, 0.25];

% Force legend markers to be square
leg.ItemTokenSize = [20, 20]; % Set size for legend markers (width, height)

% Set consistent axis limits
xlim([4.5 8.5]); % Use a fixed range for x-axis
ylim([0,max(temp.tot_scale_density)+0.0000000005]);

box on; % Adds a box around the current axes

% Adjust figure aesthetics
set(gca, 'FontSize', 36, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure
set(gcf, 'PaperPositionMode', 'auto');  % Ensures export respects figure size

fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1500, 1000];

%Turn on to save a high-resolution figure
%print(fig, 'CombinedDensity.png', '-dpng', '-r600')

hold off;
clear lab_name temp i ax leg fig
disp('Figure 1 Saved...')


%% Calculate the density of the adjusted data
% This section will calculate the adjusted lab kernel density which means
% the total density without the Intercomparison material. This is mean to
% provide a baseline comparison for the new data  to compare to the data
% which includes the intercomparison material (calculated above in the
% figure producing section.

% Remove Intercomparison data from the dataset
idxlab = ~strcmp(extlab, "Intercomparison");
CA.adjconc = CA.conc(idxlab);
CA.adjerror = CA.error(idxlab);

% Define the global x-range for the adjusted data
temp.buffer = 1 * (max(CA.adjconc) - min(CA.adjconc)); % Define the buffer
conc_range.trim_comp = linspace(min(CA.adjconc) - temp.buffer, max(CA.adjconc) + temp.buffer, 1000);

% Initialize cumulative density
temp.tot_density = zeros(1, length(conc_range.trim_comp));

% Loop through each lab to compute and stack densities
for i = 1:length(labtrim)-1
    % Add the current density to the cumulative density
    temp.stacked_density = temp.tot_density + labtrim(i).totalconc;

    % Update cumulative density
    temp.tot_density = temp.stacked_density;
end

% Store the cumulative density in the CA structure for use later
CA.adj_tot_density = temp.tot_density;

% Clear unnecessary variables from the workspace
clear temp i adjlab idxlab


%% Plot the individual lab data in a subplot
% This figure will plot the cumulative kernel density for each lab on its
% own subplot which allows for easier identification of trends within each
% lab. It also shows a violin plot (similar to a boxplot) of the data which
% shows how the mean of each lab compares to the mean of each other lab.

% Make a new figure
figure;

% Subplot positions defined manually for figure structure
subplot_positions = [
    0.09, 0.63, 0.25, 0.21;  % Top-left
    0.34, 0.63, 0.25, 0.21;  % Top-center
    0.59, 0.63, 0.25, 0.21;  % Top-right
    0.09, 0.42, 0.25, 0.21;  % Middle-left
    0.34, 0.42, 0.25, 0.21;  % Middle-center
    0.59, 0.42, 0.25, 0.21;  % Middle-right
    0.09, 0.21, 0.25, 0.21;  % Bottom-left
    0.34, 0.21, 0.25, 0.21;  % Bottom-center
    0.59, 0.21, 0.25, 0.21;  % Bottom-right
];

% Define figure subcaption lettering
temp.subcaptions = {'A. ', 'B. ', 'C. ', 'D. ', 'E. ', 'F. ', 'G. ', 'H. '};

% Predefine axes handles (gobjects returns an n x n graphics object array)
axes_handles = gobjects(length(lab) + 1, 1); % +1 makes this  3 x 3

% Run through a for loop to plot each lab after making some calculations
for i = 1:length(lab)
    % Calculate the range of values to use on the x-axis including a buffer
    lab(i).range = linspace(min(lab(i).conc) - 3*max(lab(i).error), max(lab(i).conc) + 3*max(lab(i).error), 1000);
    
    % Calculate the scaled range
    lab(i).scaledrange = lab(i).range/1e5;

    % Store the name of each lab
    lab_name = lab(i).name;

    % Store the color for each lab
    lab(i).color = color_map(lab_name);

    % Calculate the kernel density plots using a normal distribution
    lab(i).y = normpdf(lab(i).range, lab(i).conc, lab(i).error);

    % Calculate the summed probability of the densities at each site
    lab(i).totalconc = sum(lab(i).y);

    % Calculate the maximum likelihood and mean of the CRONUS-A
    % concentration for each lab
    [~, idx] = max(lab(i).totalconc);
    lab(i).maxlike = lab(i).range(idx);
    lab(i).avg = mean(lab(i).conc);
    lab(i).std = std(lab(i).conc);

    % Plotting
    axes_handles(i) = subplot('Position', subplot_positions(i, :));
    hold on;
    fill(lab(i).scaledrange, lab(i).totalconc, lab(i).color, 'LineWidth', 1,...
        'FaceAlpha', 0.6, 'EdgeColor', 'none');
    plot(lab(i).scaledrange, lab(i).y, 'Color', 'k');
    xlim([4.5, 8.5]);
    ax = gca;
    ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
    ax.XAxis.Exponent = 0;
    xticks(5:0.5:8);
    ax.XAxis.TickLabelFormat = '%.1f'; % Keep scientific notation
    ax.YTickLabel = [];
    ax.YTick = [];
    set(gca, 'LineWidth', 1, 'Box', 'on');
    text(0.05, 0.9, [temp.subcaptions{i} lab(i).name], 'Units', 'normalized', "FontSize", 32, 'FontWeight', 'normal', 'FontName', 'Helvetica');
    set(gca, 'FontSize', 32, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis

    % Add axis labels and adjust positioning conditionally
    if i == 4 %|| i == 1
        ylabel("Probability Density");
        set(gca, 'XTickLabel', [])
    elseif i == 8
        %xlabel('^{14}C Concentration (10^{5} atoms g^{-1})');
        %set(gca, 'YTickLabel', [])
    elseif i == 7
        %xlabel('^{14}C Concentration (10^{5} atoms g^{-1})'); 
        %ylabel("Probability Density");
    elseif i == 3
        xlabel('^{14}C Concentration (10^{5} atoms g^{-1})'); 
        set(gca, 'XAxisLocation', 'top');
    end
    
    % Remove axis labels for specific subplots
    if ismember(i, [2, 5, 6])
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
    end
end

% Hard code the mean and standard deviation for the University of Cologne
% based on the value from Fulop 2015
lab(8).avg = 672000;
lab(8).std = 71000;

% Super X label
annotation('textbox', [0.14, 0.07, 0.4, 0.05], ...
    'String', '^{14}C Concentration (10^{5} atoms g^{-1})', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 32, 'FontWeight', 'Normal', 'FontName', 'Helvetica', ...
    'EdgeColor', 'none');

% Super Y label
annotation('textbox', [0.01, 0.3, 0.05, 0.4], ...
    'String', 'Probability Density', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 32, 'FontWeight', 'Normal', 'FontName', 'Helvetica', ...
    'EdgeColor', 'none', ...
    'Rotation', 90);


% Add violin plots in the bottom right with adjusted position
axes_handles(end) = subplot('Position', subplot_positions(9, :));
violin(CA.matrix_data/1e5, 'xlabel', cellstr(const.labs_mod), 'facecolor', [1 1 1], 'medc', []);
ylabel("^{14}C Concentration (10^{5} atoms g^{-1})");
text(0.05, 0.9, 'I. Comparison', 'Units', 'normalized', "FontSize", 32, 'FontWeight', 'normal', 'FontName', 'Helvetica');
yticks(5:1:10);
set(gca, 'YAxisLocation', 'right');  % Move x-axis labels to the top

% Link axes for easy comparison
linkaxes(axes_handles(1:end-1), 'xy');

% Set some specific features for the figure to make it look nicer
set(gca, 'FontSize', 32, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure
set(gcf, 'PaperPositionMode', 'auto');  % Ensures export respects figure size

% Set figure size
fig = gcf; % Get current figure handle
fig.Position = [75, 75, 1500, 1200]; 

%Turn on to save a high-resolution figure
%print(fig, 'labspecific.png', '-dpng', '-r600')

hold off;
clear i idx lab_name temp ax fig subplot_positions ans axes_handles;
disp('Figure 2 Saved...')


%% Plot the individual lab data each on their own plots
% This figure will plot the cumulative kernel density for each lab on its
% own subplot which are meant to be organized in illustrator to make the
% figure in the paper, these by themselves are not as useful

% Run through a for loop to plot each lab
for i = 1:length(lab)
    % Make a new figure for each lab
    figure;

    % Store the name of each lab
    lab_name = lab(i).name;

    % Plotting
    hold on;
    fill(lab(i).scaledrange, lab(i).totalconc, lab(i).color,...
        'FaceAlpha', 0.6, 'EdgeColor', 'none');
    plot(lab(i).scaledrange, lab(i).y, 'LineWidth', 2, 'Color', 'k');
    xlim([4.5, 8.5]);
    % LDEO, Intercomparison, and ETH set to 0.00045
    % Purdue, P-CEGS, ANSTO set to 0.00015
    % TU-CEGS and Cologne set to 0.00020
    ylim([0 0.00020]);
    ax = gca;
    ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
    ax.XAxis.Exponent = 0;
    xticks(5:0.5:8);
    ax.XAxis.TickLabelFormat = '%.1f'; % Keep scientific notation
    ax.YTickLabel = [];
    ax.YTick = [];
    set(gca, 'LineWidth', 1, 'Box', 'on');
    set(gca, 'FontSize', 32, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis

    % Set some specific features for the figure to make it look nicer
    set(gca, 'FontSize', 32, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
    set(gcf, 'Color', 'w'); % White background for the figure
    set(gcf, 'PaperPositionMode', 'auto');  % Ensures export respects figure size

    % Set figure size
    fig = gcf; % Get current figure handle
    fig.Position = [75, 75, 500, 400]; 

    %Turn on to save a high-resolution figure
    filename = sprintf('LabDataFigures/labspecific_%s.png', lab_name);
    %print(fig, filename, '-dpng', '-r600');
    hold off;
end

% Reorganize the order of the labs to correspond to the individual plots
CA.matrix_data_reorg = zeros(size(CA.matrix_data));
CA.matrix_data_reorg(:,1) = CA.matrix_data(:,7);
CA.matrix_data_reorg(:,2) = CA.matrix_data(:,6);
CA.matrix_data_reorg(:,3) = CA.matrix_data(:,4);
CA.matrix_data_reorg(:,4) = CA.matrix_data(:,5);
CA.matrix_data_reorg(:,5) = CA.matrix_data(:,2);
CA.matrix_data_reorg(:,6) = CA.matrix_data(:,1);
CA.matrix_data_reorg(:,7) = CA.matrix_data(:,3);
CA.matrix_data_reorg(:,8) = CA.matrix_data(:,8);

figure;
% Create violin plots
violin(CA.matrix_data_reorg/1e5, 'facecolor', [1 1 1], 'medc', []);
%ylabel("^{14}C Concentration (10^{5} atoms g^{-1})");
%text(0.05, 0.9, 'I. Comparison', 'Units', 'normalized', "FontSize", 32, 'FontWeight', 'normal', 'FontName', 'Helvetica');
yticks(5:1:10);
set(gca, 'YAxisLocation', 'right');  % Move x-axis labels to the top

set(gca, 'LineWidth', 1, 'Box', 'on');
set(gca, 'FontSize', 32, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis

% Set some specific features for the figure to make it look nicer
set(gca, 'FontSize', 32, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure
set(gcf, 'PaperPositionMode', 'auto');  % Ensures export respects figure size

% Set figure size
fig = gcf; % Get current figure handle
fig.Position = [75, 75, 500, 400]; 

%Turn on to save a high-resolution figure
filename = sprintf('LabDataFigures/%s.png', "violin_plots");
fig = gcf; % Get current figure handle
%print(fig, filename, '-dpng', '-r600') % Turn on to save the figure

hold off;
clear i idx lab_name temp ax fig subplot_positions ans axes_handles;
disp('Lab Figures Saved...')

%% Calculate the maximum likelihood and mean of full compilation
% This section calculates the maximum likelihood and mean of the full
% composition of the data and adjusted compilation of data. This is needed
% for performing the production rate calculation because we need the
% statistically likely value.

% Calculate the statistics for the full compilation
[~, I] = max(CA.tot_density); % Index of maximum value
CA.maxlike = conc_range.full_comp(I); % Find the value corresponding to the max
CA.avgdens = mean(CA.conc);
CA.err_dens = std(CA.conc);

% Calculate the statistics for the adjusted composition (without the
% intercomparison material)
[~, Iadj] = max(CA.adj_tot_density);
CA.adjmaxlike = conc_range.trim_comp(Iadj);
CA.adjavgdens = mean(CA.adjconc);
CA.adjerr_dens = std(CA.adjconc);

clear I Iadj


%% Determine the variability of the lab data from the compilation
% This block calaculates the difference of the maximum likelihood value for
% each lab from the consensus value using two uncertainty evaluation
% methods, one that accounts for independence and one that does not

% Run error propogation method (not used, just here for evaluation)
for i=1:length(lab)
   lab(i).consensus_diff = lab(i).maxlike - CA.maxlike;
   lab(i).propogate_unc = sqrt((lab(i).std^2)+(CA.err_dens^2));
   lab(i).zscore = lab(i).consensus_diff/lab(i).propogate_unc;
end

% Leave one out mean comparison is just for evaluation, note that the
% results differ significantly from the other methods, this is due to the
% sensitivity of this test to the number of samples compared to the above
% result
for i=1:length(lab)
    lab(i).cl = ((102*CA.maxlike)-(length(lab(i).conc)*lab(i).maxlike))/(102-length(lab(i).conc));
    lab(i).cl_delta_zscore = (CA.maxlike-lab(i).cl)/CA.err_dens;
end

clear i

%% Plot the CRONUS-A Concentration through time
% This section plots the extracted value of CORNUS-A through time if the
% lab that made the measurement reports the date of the extraction. This is
% to try to evaluate if there are any trends in the CRONUS-A cncentration
% through time.


% Use the defined color_map to assign colors to the labs
temp.data_colors = cell2mat(values(color_map, cellstr(date_extracted.lab)));

% Plot the scatter plot for each unique category
figure;
hold on; % Allow multiple plots on the same figure
temp.unique_labs = unique(cellstr(date_extracted.lab)); % Get unique lab names

for i = 1:length(temp.unique_labs)
    temp.lab_idx = temp.unique_labs{i};
    temp.lab_indices = strcmp(cellstr(date_extracted.lab), temp.lab_idx); % Find indices for this lab
    if strcmp(temp.lab_idx, 'TU-CEGS')
        scatter(date_extracted.date(temp.lab_indices), date_extracted.conc(temp.lab_indices), 100, ...
            color_map(temp.lab_idx), '^', 'filled', 'DisplayName', temp.lab_idx); % Plot with the corresponding color
    else
        scatter(date_extracted.date(temp.lab_indices), date_extracted.conc(temp.lab_indices), 100, ...
            color_map(temp.lab_idx), 'o', 'filled', 'DisplayName', temp.lab_idx); % Plot with the corresponding color
    end
end

% Add legend
legend('show', 'Location', 'best'); % Automatically display all plotted data with 'DisplayName'

% Add horizontal line and annotations
yline(CA.maxlike, 'k-.', 'DisplayName', 'Max Concentration', 'FontSize', 16, 'LineWidth', 2); 

% Add text with customized font size and move it to align with the legend
text(min(date_extracted.date)+800, 700000, sprintf('CRONUS-A: %.0f (atoms g^{-1})', CA.maxlike), ...
    'FontSize', 16, 'FontWeight', 'bold', 'BackgroundColor', 'white',...
    'HorizontalAlignment','center', 'VerticalAlignment', 'middle'); % Original position

% Add labels and title
xlabel('Year');
ylabel('CRONUS-A Concentration (atoms g^{-1})');

% Add legend
legend('show', 'Location', 'best'); % Automatically display all plotted data with 'DisplayName'

% Improve formatting
set(gca, 'FontSize', 16, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure

% Optional: Adjust axis limits
xlim([min(date_extracted.date) max(date_extracted.date)]);
ylim([500000 max(date_extracted.conc) * 1.1]); % Adjust limits for better visibility

% Set figure size
fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1200, 800]; % Example: 1200x800 pixels

%Turn on to save a high-resolution figure
%print(fig, 'cronusa_date.png', '-dpng', '-r600')

hold off;
disp('Figure 3 Saved...')

clear fig i lab_indices lab_idx


%% Assess the statistical similarity of the lab measurements
% In order to robustly determine the similarity of the labs extraction
% values I will perform a one-way Anova to test if they are statistically
% similar or not.

% Preform the Anova on the raw data to see if it is different
[anova.p, anova.tbl, anova.stats] = anova1(CA.conc_raw, CA.extlab_raw);

% Compare each lab individually, this shows if eah lab is statistically
% different from the others based on the results of the Anova.
figure;
anova.restuls = multcompare(anova.stats);


% Display a progress report
disp('ANOVA completed...')


%% Establish Production Rate Vector
% Establish a range of SLHL production rate values to try based on
% an estimated range from the literature.

% The minimum value is set to 7 to accomodate a low value in the
% uncertainty of spallation production after accounting for muon
% production. The general range of values from the literature for
% spallation is around 9-16 atoms/g/yr.

% To extend the range of values or decrease the spacing (default is 0.001
% atoms/g/yr) modify the min, max and interval variables below.

P14.min = 7; % minimum estimate
P14.max = 20; % maximum estimate
P14.intervals = (P14.max-P14.min)*1300; % number of interval values to check
P14.range = linspace(P14.min, P14.max, P14.intervals); % vector of possible production rates
 

%% Stone Scaling Production Rate Estimate
% This code block runs the StPRodRate function which performs a production
% rate calibration of an input site (in this case the CRONUS-A site) and
% calibrates the production rate for each lab where there is data using the
% Stone, 2000 scaling scheme to scale the estimate to SLHL and the ERA40
% atmospheric model. 

% Write a progress report
disp('Calculating the Stone Scaling Production Rate...')

% Run the StProdRate function which calculates the production rate of
% CRONUS-A using the Stone scaling framework
[St, lab] = StProdRate(CA, 'CRONUSA', lab, P14, const);


%% Monte-Carlo estimate of CRONUS-A distribution
% This section performs a monte-carlo simulation of the CRONUS-A
% distribution to create a 1x10000 vector of likely CRONUS-A concentrations
% that can be used to determine more robustly the uncertainty on the
% production rate. The matrix is input as a matrix into the production rate
% calculation code in order to produce an equal size matrix of production
% rate values

% Write a progress report
disp('Beginning the Monte-Carlo Simulation...')

% Initialize the for loop
for i=1:10000
    % Make sure we actually get 10000 points 
    while true
        % Draw a random value with uniform probability over the full range of
        % the CRONUS-A data
        temp.x_rand(i) = min(conc_range.full_comp) + (max(conc_range.full_comp)-min(conc_range.full_comp))*rand;

        % Calculate the PDF probability for that value using a linear
        % interpolation
        temp.y_rand(i) = interp1(conc_range.full_comp, CA.tot_density, temp.x_rand(i), 'linear', 0);

        % Use a conditional to accept or not accept the value based on the pdf
        % and store the value in an array
        if rand*max(CA.tot_density) < temp.y_rand(i)
            CA.MC(i) = temp.x_rand(i);
            break; % Exit the while loop
        end
    end
end

clear temp i


%% Uncertainty using Monte-Carlo Simulation Results (St Scaling)
% This section takes the Monte Carlo results from the previous code section
% which represent 10000 likely CRONUS-A samples and calculates the
% production rate for CRONUS-A assuming that every one of these was the
% correct value. This is the least inclusive of the multiple different
% methods I used to evaluate the uncertainty.

 % Add a progress report
 disp('Calculating the Monte Carlo Production Rates using St Scaling...')

 % For each value in the Monte Carlo Simulation calculate the SLHL
 % production rate using the P14_map established earlier
 for i=1:length(CA.MC)
    [~,temp.idx] = min(abs(St.Nsat-CA.MC(i)));
    temp.key(i) = St.Nsat(temp.idx);
    St.MC_Prod_rate(i) = St.P14_map(temp.key(i));
 end

 % Sort the values in ascending order
 St.MC_Prod_rate = sort(St.MC_Prod_rate, 'ascend');

 % Calculate the kernel density of the Monte Carlo production rates
 St.MC_kde_dens = kde(St.MC_Prod_rate, "NumPoints", 10000);

 % Get the index of the maximum kernel density value
 [~,temp.I] = max(St.MC_kde_dens);

 % Calculate the maximum likelihood and standard deviation of the Monte
 % Carlo Production Rates
 St.MC_maxlike_Prod_rate =  St.MC_Prod_rate(temp.I);
 St.MC_stddev_Prod_rate = std(St.MC_Prod_rate);

 clear temp i

%% Uncertainty Analysis using Parameter Ranges (St scaling)
% In this section an alternative method to evaluate the uncertainty on the
% production rate is used by altering some of the input values to the
% CRONUS-A produciton rate claibration in order to calculate the lower
% possible range based on the parameter uncertainty. Assumes the same muon
% production rate as the main calculation, this only varies the other
% parameters within their uncertainty so is the second most robust estimate
% of uncertainty.

% Write a progress report
disp('Calculating the Stone Scaling Production Rate Uncertainty Bounds...')

% Genereate a new structure of data for the two ranges of uncertainty
% analysis
CA.CAmin.lat = CA.lat; CA.CAmin.long = CA.long; CA.CAmin.z = CA.z+CA.z_error; CA.CAmin.maxlike = CA.maxlike-CA.err_dens; CA.CAmin.avgdens = CA.avgdens-CA.err_dens;
CA.CAmax.lat = CA.lat; CA.CAmax.long = CA.long; CA.CAmax.z = CA.z-CA.z_error; CA.CAmax.maxlike = CA.maxlike+CA.err_dens; CA.CAmax.avgdens = CA.avgdens+CA.err_dens;

% Run the StProdRate function which calculates the production rate of
% CRONUS-A using the Stone scaling framework
[St.Minunc, lab] = StProdRate(CA.CAmin, 'CRONUSA', lab, P14, const.constmin);
[St.Maxunc, lab] = StProdRate(CA.CAmax, 'CRONUSA', lab, P14, const.constmin);


%% Full-variability Monte-carlo (St scaling)
% In this section a different Monte-carlo procedure is performed which
% varies each of the three variables (decay constant, elevation, concentration)
% that has the potential to impact the
% production rate, then determines the relevant production rate to generate
% a full distribution that includes the full range of uncertainty.

%NOTE: This section takes a while (10 minutes for parallel processing with
% 14 processesors so unless you are planning to evalaute a more robust uncertainty 
% than can be accomplished in the lines above, leave this commented out for 
% the sake of time. If you do want to run it then you will need to have the parallel processing package
% installed, and will need to change the M value (in the parfor loop after
% the variable range) to establish the maximum number of processors that can
% be used based on your available compute power on your own machine

% COMMENT OR UNCOMMENT FROM THE LINE BELOW
% Write a progress report
% disp('Calculating the Stone Scaling Production Rate Range using a Monte-Carlo...')
% 
% % Set up the distribution before running the loop
% % Normalize the distribution and calculate the cdf 
% cdf = cumsum(CA.tot_density/sum(CA.tot_density));
% 
% MC_local = zeros(1,10000);
% 
% % It is possible to run the same simulation with a smaller interval which
% % cuts the computation time down by 2 orders of magnitude (using this
% % adjusted range gives the exact same value but takes less than10 minutes).
% P14_adj.range = [7:0.02:20];
% 
% % Evaluate the time needed to run this
% tic;
% 
% parfor (i=1:10000,14)
%     %Initialize variables
%     temp = struct();
%     ctemp = struct();
% 
%     % Set constants
%     temp.lat = CA.lat;
%     temp.long = CA.long;
%     ctemp.Fsp = 1;
% 
%     % Draw a random elevation
%     temp.z = random('Normal', CA.z, CA.z_error);
% 
%     % Draw a random decay constant
%     temp.halflife = random('Normal',5700,30);
%     ctemp.lambda = log(2)./temp.halflife;
% 
%     % Draw a random CRONUS-A concentration
%     r = rand(1,1);
%     temp.conc = interp1(cdf, conc_range.full_comp, r, 'linear','extrap');
% 
%     % Draw a random muon reference production rate assuming that the mean
%     % is the average of an empirically derived and experimentally
%     % determined value and the standard deviation between the two
%     temp.muon_ref = random('Normal', 3.43, 0.5);
%     temp.muonPR = temp.muon_ref*exp((1013.25-const.atmp)/const.L)
% 
%     MC_local(i) = StProdRateUnc(temp, P14_adj, ctemp);
% end
% 
% % Save the local allocation to the full dataset
% St.MC_fullvariability = MC_local;
% 
% %Print the time needed to run this section
% elapsedTime = toc;
% fprintf('Elapsed time for this section: %.4f seconds\n', elapsedTime);
% 
% clear temp ctemp cdf r i MC_local elapsedTime
% COMMENT OR UNCOMMENT TO THE LINE ABOVE


% If you choose to run the above section then comment out these lines of code
% below so as not to overwrite your data
load('MC_fullvariability_St_muons.mat')
St.MC_fullvariability = MC_St_muons;
clear MC_St_muons

% Calculate the uncertainty based on the full variability Monte Carlo
% Experiment
St.MC_fullvariability_SD = std(St.MC_fullvariability);
St.MC_fullvariability_avg = mean(St.MC_fullvariability);

%% Plot full variability Monte-carlo for St scaling
% Plot the data from above to visualize the distribution of production
% rates that is generated
figure;

histogram(St.MC_fullvariability, 'Normalization','probability');
xlabel('Production Rate', 'Fontsize', 24, 'FontName', 'Helvetica', 'Fontweight', 'Normal');
ylabel('Probability', 'Fontsize', 24, 'FontName', 'Helvetica', 'Fontweight', 'Normal');

% Improve formatting
set(gca, 'FontSize', 24, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure


%% LSDn Scaling Section
% This section calculates the time dependent production rate using the LSDn 
% Scaling (Lifton et al., 2014). The code here is from the CD14C code published
% in Koester and Lifton, 2023 which calculates the compositionally
% dependent production rate using a P14 of 13.5 for quartz from
% spallation. This version is modified to take a range of possible
% production values (corresponding to the same range used in the St scaling
% method) to output a range of time dependent production rates that can be
% used to determine the best fit production rate.

% Display a progress report
disp('Calculating the LSDn scaled production rate...');
    
LSD = CD14C_CRONUScalib('LSDn_inputs.txt', 1);
% The output here "LSD" is a structure that stores the time and
% compostitionally dependent production rate at specific time intervals
% given by the vector "tv" defining the time period of each production rate
% estimate.


%% Calculate the saturation concentration using LSDn scaling
% First find the indices of the time vector less than 50,000 years. We expect any
% continuously exposed material that is not undergoing erosion to be satured
% with respect to in-situ 14C after around 30,000 years, so 50,000 years is
% a conservative estimate. This step reduces processing time.

% Add a progress report
disp('Calculating the CRONUS-A Saturation Concentration...')

temp.idx = find(LSD.tv<500000); % Find the indices of the time <50,000 years
temp.sattv = LSD.tv(temp.idx); % Get the actual times for those indices

% Next run through the output of the CD14C code to find the production rate
% for each timestep corresponding to the trimmed time vector.
for i=1:size(LSD.P14_CD, 1)
    temp.satP14(i,:) = LSD.P14_CD(i,temp.idx);
end

% Flip both vectors so that the accumulation starts 50 ka instead of in the
% present.
LSD.tv = flip(temp.sattv);
LSD.P14v = fliplr(temp.satP14);


% Calculate the muon production rate to subtract from the total
% reference production rate. Assumes average production rate of Balco
% 2017 (3.07 atoms/g/yr) and Heisinger 2002 (3.78 atoms/g/yr) for total
% muon SLHL surface production (uses 3.43 atoms/g/yr)
muon_PR = 3.43.*exp((1013.25-ERA40atm(CA.lat, CA.long, CA.z))/const.L);



% Eqution 4.1 of Dunai, 2010 assuming that the initial inventory C_inh is
% zero calculates the saturation calculation for cosmogenic nuclides
for b=1:size(LSD.P14v, 1)
    for a=1:length(LSD.tv)
        LSD.C14sum(b,a) = ((1 - exp(-1*LSD.tv(a).*const.lambda)).*(LSD.P14v(b,a)+muon_PR))./const.lambda;
    end
end


% Store the maximum value in a structured array for each access later
LSD.max = max((LSD.C14sum)');

clear a b i temp;


%% Plot the time-dependent 14C accumulation
% This plot shows the accumulation of cosmogenic nuclides up to the
% saturation concentration using the equation 4.1 from Dunai, 2010. The
% plot is meant to demonstrate the relationship.

% Plot the log log relationship because the features are on the log scale
figure;
loglog(LSD.tv, LSD.C14sum(8293,:), 'LineWidth', 3,...
       'DisplayName', '^{14}C concentration curve', 'Color', 'k');

% Improve grid and formatting
grid on; % Turn on the grid
set(gca, 'MinorGridLineStyle', '-'); % Set gridlines to solid

% Customize the grid lines to make them less bright
ax = gca; % Get the current axes
ax.GridAlpha = 0.2; % Set grid transparency (lower value for less brightness)
ax.GridColor = [0.5, 0.5, 0.5]; % Set grid color (lighter gray)
ax.MinorGridAlpha = 0.05;
ax.MinorGridColor = [0.5, 0.5, 0.5];

% Customize axes labels and make them larger
xlabel('Exposure duration (years)');
ylabel('In-situ ^{14}C concentration (atoms g^{-1})');

% Add a horizontal line at the saturation concentration
saturation_concentration = round(LSD.max(8293)); % Get the saturation concentration as an integer
yline(saturation_concentration, '-.r', 'LineWidth', 3, ...
    'DisplayName', 'Saturation Concentration'); % Dashed red line for visibility

% Add text with customized font size and move it to align with the legend
text(40, 500000, sprintf('Saturation: %d atoms g^{-1}', saturation_concentration), ...
    'FontSize', 24, 'FontWeight', 'normal', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1.5); % Original position

% Add a legend with both the curve and horizontal line
legend('show', 'Location', 'best', 'FontSize', 24); % Align legend near the text

% Adjust figure aesthetics
set(gca, 'FontSize', 24, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure

% Set figure size
fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1200, 800]; % Example: 1200x800 pixels

%Turn on to save a high-resolution figure
%print(fig, 'age_vs_acc_to_sat.png', '-dpng', '-r600')

clear saturation_concentration ax 

disp('Figure 4 Saved...')


%% Estimate the LSDn Reference Production Rate
% This section will use the saturation concentrations that were calculated
% above to establish a dictionary of saturation values and production rates
% that determine those values, then it will calculate the best fit of the
% production rate based on the saturation values.

% Add a progress report
disp('Calculating the LSDn reference production rate...')

% Create a dictionary of the saturation values and the reference production
% rates that result in those saturation values
LSD.P14_map = dictionary(LSD.max, P14.range);

% Find the saturation value of the elevation scaled value that is closest 
% to the CRONUS-A measurement from all of the compiled data and then 
% determines what the SLHL scaled value is.
[~,temp.avgidxLSD] = min(abs(LSD.max-CA.avgdens));
[~,temp.maxidxLSD] = min(abs(LSD.max-CA.maxlike));
temp.avgkeyLSD = LSD.max(temp.avgidxLSD);
temp.maxkeyLSD = LSD.max(temp.maxidxLSD);
LSD.P14SLHLmax = LSD.P14_map(temp.maxkeyLSD);
LSD.P14SLHLavg = LSD.P14_map(temp.avgkeyLSD);

% Perform the same calculation to determine the SLHL scaled production rate
% from the saturation value for each lab measurement and store those values
% in the relevant laboratory structure
for i=1:length(lab)
    [~,lab(i).avgidx] = min(abs(LSD.max-lab(i).avg));
    [~,lab(i).maxidx] = min(abs(LSD.max-lab(i).maxlike));
    lab(i).avgkey = LSD.max(lab(i).avgidx);
    lab(i).maxkey = LSD.max(lab(i).maxidx);
    lab(i).LSDP14SLHLmax = LSD.P14_map(lab(i).maxkey);
    lab(i).LSDP14SLHLavg = LSD.P14_map(lab(i).avgkey);
end

clear i temp;


%% Use the Monte Carlo Simulation to Calculate LSDn Uncertainty
% Similar to the procedure performed above for the St scaling method I am
% going to determine the uncertainty of the LSDn production rate using a
% monte carlo simulation of the likely CRONUS-A values. I implement the
% same monte carlo simulation as above for the sake of consistency and
% comparison between the production rates.

 % Add a progress report
 disp('Calculating the Monte Carlo Production Rates using LSDn Scaling...')


 % For each value in the Monte Carlo Simulation calculate the SLHL
 % production rate using the P14_map established earlier
 for i=1:length(CA.MC)
    [~,temp.idx] = min(abs(LSD.max-CA.MC(i)));
    temp.key(i) = LSD.max(temp.idx);
    LSD.MC_Prod_rate(i) = LSD.P14_map(temp.key(i));

 end

 % Sort the values in ascending order
 LSD.MC_Prod_rate = sort(LSD.MC_Prod_rate, 'ascend');

 % Calculate the kernel density of the Monte Carlo production rates
 LSD.MC_kde_dens = kde(LSD.MC_Prod_rate, "NumPoints", 10000);

 % Get the index of the maximum kernel density value
 [~,temp.I] = max(LSD.MC_kde_dens);

 % Calculate the maximum likelihood and standard deviation of the Monte
 % Carlo Production Rates
 LSD.MC_maxlike_Prod_rate =  LSD.MC_Prod_rate(temp.I);
 LSD.MC_stddev_Prod_rate = std(LSD.MC_Prod_rate);
 
 clear temp i


 %% Full-variability Monte-carlo (LSDn scaling)
% In this section a different Monte-carlo procedure is performed which
% varies each of the three variables (decay constant, elevation, concentration)
% that has the potential to impact the
% production rate, then determines the relevant production rate to generate
% a full distribution that includes the full range of uncertainty.

%NOTE: This section takes an extremely long time to run (3 hours or more for parallel processing
% so unless you are planning to evalaute a more robust uncertainty than can be accomplished in
% the lines above, leave this commented out for the sake of time. If you do
% want to run it then you will need to have the parallel processing package
% installed, and will need to change the M value (in the parfor loop after
% the variable range) to establish the maximum number of processors that can
% be used based on your available compute power on your own machine

% COMMENT OR UNCOMMENT FROM THE LINE BELOW
% % Write a progress report
% disp('Calculating the LSDn Scaling Production Rate Range using a Monte-Carlo...')
% 
% % Set up the distribution before running the loop
% % Normalize the distribution and calculate the cdf 
% cdf = cumsum(CA.tot_density/sum(CA.tot_density));
% 
% MC_local = zeros(1,10000);
% P14.adjrange = [7:0.02:20];
% 
% % Evaluate the time needed to run this
% tic;
% 
% parfor (i=1:10000,14)
%     %Initialize variables
%     temp = struct();
%     itemp = struct();
% 
%     % Set constants
%     temp.name = 'CRONUS-A';
%     temp.lat = CA.lat;
%     temp.long = CA.long;
%     temp.xrf = [100 zeros(1,10)];
% 
%     % Draw a random elevation
%     temp.z = random('Normal', CA.z, CA.z_error);
% 
%     % Draw a random decay constant
%     temp.halflife = random('Normal',5700,30);
%     temp.lambda = log(2)./temp.halflife;
% 
%     % Draw a random CRONUS-A concentration
%     r = rand(1,1);
%     itemp.conc = interp1(cdf, conc_range.full_comp, r, 'linear','extrap');
% 
%     % Draw a random muon reference production rate assuming that the mean
%     % is the average of an empirically derived and experimentally
%     % determined value and the standard deviation between the two
%     temp.muon_ref = random('Normal', 3.43, 0.5);
%     temp.muonPR = temp.muon_ref*exp((1013.25-const.atmp)/const.L)
% 
% 
%     MC_local(i) = LSDProdRateUnc(temp, 3, itemp, P14.adjrange)
% end
% 
% % Save the local allocation to the full dataset
% LSD.MC_fullvariability = MC_local;
% 
% %Print the time needed to run this section
% elapsedTime = toc;
% fprintf('Elapsed time for this section: %.4f seconds\n', elapsedTime);
% 
% clear temp itemp cdf r i MC_local elapsedTime
% COMMENT OR UNCOMMENT TO THE LINE ABOVE


% If you choose to run the above section then comment out these lines of code
% below so as not to overwrite your data
load('MC_fullvariability_LSDn_muons.mat')
LSD.MC_fullvariability = MC_LSDn_muons;
clear MC_fullvariability

% Calculate the uncertainty based on the full variability Monte Carlo
% Experiment
LSD.MC_fullvariability_SD = std(LSD.MC_fullvariability);
LSD.MC_fullvariability_avg = mean(LSD.MC_fullvariability);

%% Plot full variability Monte-carlo using LSDn scaling
% Plot the data from above to visualize the distribution of production
% rates that is generated
figure;

histogram(LSD.MC_fullvariability, 'Normalization','probability');
xlabel('Production Rate', 'Fontsize', 24, 'FontName', 'Helvetica', 'Fontweight', 'Normal');
ylabel('Probability', 'Fontsize', 24, 'FontName', 'Helvetica', 'Fontweight', 'Normal');

% Improve formatting
set(gca, 'FontSize', 24, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure



%% Calculate the Elevation Scaled Saturation Curves
% This section calculates the saturation curve specific for the CRONUS-A
% site at 31 different elevations to show how the saturation value of a
% site changes based on elevation. 

% Add a progress report
disp('Calculating the Elevation Scaled Saturation Curves...')

z_range=0:100:3000; %elevation range to calculate the values

% For loop that calculates the scaled production rate at elevation based 
% on estimated SLHL production rates calculated using different statistical
% likelihoods and the St scaling framework. It also calculates the muon
% production rate at elevation which is needed for the steps below
for i=1:length(z_range)
    St.P14avg(i,1) = St.P14SLHLavg .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp);
    St.P14max(i,1) = St.P14SLHLmax .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp);
    muons.atmp_h(i) = ERA40atm(CA.lat, CA.long, z_range(i));
    muons.elev_PR(i) = 3.43.*exp((1013.25-muons.atmp_h(i))/const.L);

    for j=1:length(lab)
        lab(j).P14avg(i,1) = lab(j).StP14SLHLavg .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp);
        lab(j).P14max(i,1) = lab(j).StP14SLHLmax .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp);
    end
end

muons.elev_PR = muons.elev_PR';

% Calculate the saturation curve concentration for different elevations
% scaled based on the production rate at those elevations.
St.Nsatavg = (St.P14avg+muons.elev_PR)./(const.lambda);
St.Nsatmax = (St.P14max+muons.elev_PR)./(const.lambda);

% Calculate the saturation curve for each lab using the same method as
% above.
for i=1:length(lab)
    lab(i).Nsatavg = (lab(i).P14avg+muons.elev_PR)./(const.lambda);
    lab(i).Nsatmax = (lab(i).P14max+muons.elev_PR)./(const.lambda);
end

% Generate a cell structure to calculate the production rate at different
% elevations for the CRONUS-A concentration. 
for i=1:length(z_range)
    LSD.satcurvein{i,1} = 'CRONUS-A';
    LSD.satcurvein{i,2} = CA.lat;
    LSD.satcurvein{i,3} = CA.long;
    LSD.satcurvein{i,4} = z_range(i);
    LSD.satcurvein{i,5} = 100;
    for ii=1:10
        LSD.satcurvein{i,ii+5} = 0;
    end
end


% Run the function to calculate the saturation concentration for each
% elevation using the LSDn scaling framework.
LSD.satcurveout = CD14C_CRONUScalib(LSD.satcurvein, 2, LSD.P14SLHLmax);

% Perform calculations to find the value of time vector less than 50 kyr.
temp.idx = find(LSD.satcurveout.tv<50000);
temp.sattv = LSD.satcurveout.tv(temp.idx); %Find the values from the indices

% Find the production rates trimmed to 50 kyr.
for i=1:size(LSD.satcurveout.P14_CD, 1)
    temp.satcurveP14(i,:) = LSD.satcurveout.P14_CD(i,temp.idx);
end

%Flip both vectors so that the accumulation starts at 50 ka not the
%present.
LSD.satcurveout.tvt = flip(temp.sattv);
LSD.satcurveout.P14v = fliplr(temp.satcurveP14);

% Calculate the saturation concentration for each elevation of the CRONUS-A
% site using the same structure as above to calculate the reference
% production rate.
for b=1:size(LSD.satcurveout.P14v, 1)
    for a=1:length(LSD.satcurveout.tvt)
        LSD.satcurveout.C14sum(b,a) = ((1 - exp(-1*LSD.satcurveout.tvt(a).*const.lambda)).*(LSD.satcurveout.P14v(b,a)+muons.elev_PR(b)))./const.lambda;
    end
end

% Find the maximum value for each accumulation that corresponds to the
% saturation value.
LSD.satcurve = max((LSD.satcurveout.C14sum)');

% Perform the same procedure for the uncertainty bounds
LSD.satcurveout_min = CD14C_CRONUScalib(LSD.satcurvein, 2, LSD.P14SLHLmax-LSD.MC_fullvariability_SD);
LSD.satcurveout_max = CD14C_CRONUScalib(LSD.satcurvein, 2, LSD.P14SLHLmax+LSD.MC_fullvariability_SD);


% Minimum
% Perform calculations to find the value of time vector less than 50 kyr.
temp.idx_min = find(LSD.satcurveout_min.tv<50000);
temp.sattv_min = LSD.satcurveout_min.tv(temp.idx_min); %Find the values from the indices

% Find the production rates trimmed to 50 kyr.
for i=1:size(LSD.satcurveout_min.P14_CD, 1)
    temp.satcurveP14_min(i,:) = LSD.satcurveout_min.P14_CD(i,temp.idx_min);
end

%Flip both vectors so that the accumulation starts at 50 ka not the
%present.
LSD.satcurveout_min.tvt = flip(temp.sattv_min);
LSD.satcurveout_min.P14v = fliplr(temp.satcurveP14_min);

% Calculate the saturation concentration for each elevation of the CRONUS-A
% site using the same structure as above to calculate the reference
% production rate.
for b=1:size(LSD.satcurveout_min.P14v, 1)
    for a=1:length(LSD.satcurveout_min.tvt)
        LSD.satcurveout_min.C14sum(b,a) = ((1 - exp(-1*LSD.satcurveout_min.tvt(a).*const.lambda)).*(LSD.satcurveout_min.P14v(b,a)+muons.elev_PR(b)))./const.lambda;
    end
end

% Find the maximum value for each accumulation that corresponds to the
% saturation value.
LSD.satcurve_min = max((LSD.satcurveout_min.C14sum)');

% Maximum
% Perform calculations to find the value of time vector less than 50 kyr.
temp.idx_max = find(LSD.satcurveout_max.tv<50000);
temp.sattv_max = LSD.satcurveout_max.tv(temp.idx_max); %Find the values from the indices

% Find the production rates trimmed to 50 kyr.
for i=1:size(LSD.satcurveout_max.P14_CD, 1)
    temp.satcurveP14_max(i,:) = LSD.satcurveout_max.P14_CD(i,temp.idx_max);
end

%Flip both vectors so that the accumulation starts at 50 ka not the
%present.
LSD.satcurveout_max.tvt = flip(temp.sattv_max);
LSD.satcurveout_max.P14v = fliplr(temp.satcurveP14_max);

% Calculate the saturation concentration for each elevation of the CRONUS-A
% site using the same structure as above to calculate the reference
% production rate.
for b=1:size(LSD.satcurveout_max.P14v, 1)
    for a=1:length(LSD.satcurveout_max.tvt)
        LSD.satcurveout_max.C14sum(b,a) = ((1 - exp(-1*LSD.satcurveout_max.tvt(a).*const.lambda)).*(LSD.satcurveout_max.P14v(b,a)+muons.elev_PR(b)))./const.lambda;
    end
end

% Find the maximum value for each accumulation that corresponds to the
% saturation value.
LSD.satcurve_max = max((LSD.satcurveout_max.C14sum)');


%% Calculate the accumulation isochrons using the LSDn method
% This section will use the satuartion output predicted by the ebst fit
% production rate and the LSDn scaling method to determine the elevation
% dependent isochrons for each sample which can be plotted below.

% Calculate the accumulation expected in isochrons that have been exposed
% for up to 20 ka
iso.t = [2000:2000:9000, 10000:5000:20000];  % Store the time intervals that we want to calculate isochrons

for k=1:length(z_range) % For each elevation interval calculate the isochrons
    iso.isobars(k,:) = interp1(LSD.satcurveout.tvt, LSD.satcurveout.C14sum(k,:), iso.t); %Interpolate the 14C between two times
end

% Clear unnecessary variables
clear k

%% Plot the saturation curves and all Antarctic 14C data from ICE-D
% This section plots the saturation cuves for the CRONUS-A site as well as
% the CRONUS-A and other saturation site 14C concentrations. All available
% CRONUS-A data is plotted as well.

% Save the data for the other saturated sites for plotting
Ant_Sat.conc = [183030; 968970; 160050; 974370; 1177930; 1038010; 1306000];
Ant_Sat.lat = [-70.86; -77.75; -70.82; -77.75; -73.44; -73.39; -77.86811];
Ant_Sat.long = [68.13; 160.8; 68.17; 160.8; 61.9; 61.72; 159.53455];
Ant_Sat.elev = [225; 2160; 100; 2020; 2538; 2137; 2696];
Ant_Sat.unc = [8420; 15770; 12860; 19180; 19490; 20640; 16830];
Ant_Sat.names = {"98-PCM-010-SRDK"; "WBC-UVP"; "98-PCM-002-BVLK"; "WBC-2020"; "98-PCM-105-MNZ"; "98-PCM-067-MNZ"; ""};

load all_Antarctic_14C.txt %text file with all in situ 14C measurements in ICE-D in Antarctica (as of 3 May 2024)
ant.z=all_Antarctic_14C(:,1); 
ant.conc=all_Antarctic_14C(:,2); 
ant.err=all_Antarctic_14C(:,3);

% Initialize figure
figure;
hold on % Plotting multiple things

% Plot the saturation curve uncertainty
LSD.curveunc = fill([LSD.satcurve_min/1e5 fliplr(LSD.satcurve_max/1e5)], [z_range fliplr(z_range)], [0.529, 0.808, 0.980], 'FaceAlpha', 0.44, 'LineStyle', 'none');

% Plot isochrons
for i=1:length(iso.t)
  iso.plt = plot(iso.isobars/1e5,z_range,':k', LineWidth=1.5, DisplayName='LSDn isochrons'); % Uses the LSDn method which is why they look off from the St scaling
end

% Add text for the isobars 
text(3.9, 3060, '2ka', 'FontSize', 24, 'FontWeight','normal','FontName','Helvetica');
text(6.9, 3060, '4ka', 'FontSize', 24, 'FontWeight','normal','FontName','Helvetica');
text(10.2, 3060, '6ka', 'FontSize', 24, 'FontWeight','normal','FontName','Helvetica');
text(12.7, 3060, '8ka', 'FontSize', 24, 'FontWeight','normal','FontName','Helvetica');
text(14.2, 3060, '10ka', 'FontSize', 24, 'FontWeight','normal','FontName','Helvetica');
text(16.5, 3060, '15ka', 'FontSize', 24, 'FontWeight','normal','FontName','Helvetica');
text(18.0, 3060, '20ka', 'FontSize', 24, 'FontWeight','normal','FontName','Helvetica');

% Plot the saturation curves and adjust the colors
St.curvemax = plot(St.Nsatmax/1e5, z_range, DisplayName='St saturation curve');
LSD.curvemax = plot(LSD.satcurve/1e5, z_range, DisplayName='LSDn saturation curve');
set(St.curvemax, 'LineWidth', 3, 'color',[0.937, 0.502, 0.502]);
set(LSD.curvemax, 'LineWidth', 4, 'color',[0.529, 0.808, 0.980]);

% Plot the CRONUS-A value with errorbars
CA.errbar = errorbar(CA.maxlike/1e5, CA.z, CA.err_dens/1e5, 'horizontal', 'LineWidth', 2, 'Color', 'k');
CA.plt = plot(CA.maxlike/1e5,CA.z,'o', MarkerSize= 22, markerfacecolor =[0.882, 0.745, 0.416], MarkerEdgeColor='k', DisplayName='CRONUS-A maximum likelihood');

% Plot all in situ 14C measurements from Antarctica
ant.curve = scatter(ant.conc/1e5, ant.z, 75, 'k', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none', DisplayName='Antarctic data');

% Plot the other Antarctic Saturated Surfaces
Ant_Sat.errbar = errorbar(Ant_Sat.conc/1e5, Ant_Sat.elev, Ant_Sat.unc/1e5, 'horizontal', 'LineWidth', 2, 'Color', 'k', 'LineStyle','none');
Ant_Sat.plt = plot(Ant_Sat.conc/1e5, Ant_Sat.elev, 'o', MarkerSize= 16, markerfacecolor = [0.251, 0.690, 0.651], MarkerEdgeColor='k', DisplayName='Other saturated surfaces');

% Labels
xlabel('^{14}C concentration (10^{5} atoms g^{-1})');
ylabel('Elevation (m asl)');
xlim([0 20]);
ylim([0 3000]);
legend([LSD.curvemax(1), St.curvemax(1), iso.plt(1), CA.plt, Ant_Sat.plt, ant.curve], 'Box', 'off');
legend('Position', [0.62, 0.3, 0.1, 0.1]);
set(gca, 'XTick', 0:2.5:20);

hold off;

% Adjust figure aesthetics
set(gca, 'FontSize', 30, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure
box on;

fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1200, 800];

%Turn on to save a high-resolution figure
%print(fig, 'saturationcurve.png', '-dpng', '-r600')

disp("Figure 5 saved...")

clear fig i

%% Calculate the saturation values for each saturated test site in Antarctica
% This section uses the calibrated reference production rate to determine
% the saturation value of the saturated sites in Antarctica in order to
% perform a chi-squared test to see how well this value fits to other data.

% Add a progress report
disp("Fitting the test saturation samples...")

for i=1:length(Ant_Sat.elev)
    Ant_Sat.LSDinput{i,1} = Ant_Sat.names(i);
    Ant_Sat.LSDinput{i,2} = Ant_Sat.lat(i);
    Ant_Sat.LSDinput{i,3} = Ant_Sat.long(i);
    Ant_Sat.LSDinput{i,4} = Ant_Sat.elev(i);
    Ant_Sat.LSDinput{i,5} = 100;
    for ii=1:10
        Ant_Sat.LSDinput{i,ii+5} = 0;
    end
end

Ant_Sat.LSD_output = CD14C_CRONUScalib(Ant_Sat.LSDinput, 2, LSD.P14SLHLmax);

temp.idx = find(Ant_Sat.LSD_output.tv<500000);
temp.sattv = Ant_Sat.LSD_output.tv(temp.idx); %Find the values

for i=1:size(Ant_Sat.LSD_output.P14_CD, 1)
    Ant_chisquare(i,:) = Ant_Sat.LSD_output.P14_CD(i,temp.idx); %Find the production rates for the time
end

%Flip both vectors so that the accumulation starts at 50 ka instead of in the
%present.
Ant_Sat.LSD_output.tvt = flip(temp.sattv);
Ant_Sat.LSD_output.P14v = fliplr(Ant_chisquare);

for i=1:length(Ant_Sat.elev)
    Ant_Sat.atmp_h(i) = ERA40atm(Ant_Sat.lat(i), Ant_Sat.long(i), Ant_Sat.elev(i));
    Ant_Sat.muon_PR(i) = 3.43.*exp((1013.25-Ant_Sat.atmp_h(i))/const.L);
end

for b=1:size(Ant_Sat.LSD_output.P14v, 1)
    for a=1:length(Ant_Sat.LSD_output.tvt)
        LSDAnt_chisquare.C14sum(b,a) = ((1 - exp(-1*Ant_Sat.LSD_output.tvt(a).*const.lambda)).*(Ant_Sat.LSD_output.P14v(b,a)+Ant_Sat.muon_PR(b)))./const.lambda;
    end
end

LSDAnt_chisquare.sat = max((LSDAnt_chisquare.C14sum)')';

clear a b i ii temp

%% Calculate chi squared values
% Calculate the chi-squared value for the saturated sites comparing the
% predicted to the measured concentrations to determine fit.

% Calculate the chi-squared goodness of fit statistic
chi2.chi_squared = sum(((Ant_Sat.conc - LSDAnt_chisquare.sat).^2)./Ant_Sat.unc.^2);
chi2.dof = 6;

% Calculate the p-value
chi2.p_val = 1 - chi2cdf(chi2.chi_squared, chi2.dof);

% Print the results 
fprintf('The chi-squared value is %.4f and the p-value is %.4f\n', chi2.chi_squared, chi2.p_val);

% This gives a p-value of 0 and a chi-squared value of 130.2 which is a
% bad fit, but there are two points that plot clearly below saturation,
% which tells me that maybe they are causing this because the other 4
% samples are at saturation. 


%% Calculate chi squared values adjusted for the possibly non-saturated samples
% Calculate the chi-squared value for the saturated sites comparing the
% predicted to the measured concentrations to determine fit.

Ant_Sat.conc_exc = Ant_Sat.conc([1,3,4,6]);
Ant_Sat.unc_exc = Ant_Sat.unc([1,3,4,6]);

LSDAnt_chisquare.sat_exc = LSDAnt_chisquare.sat([1,3,4,6]);

% Calculate the chi-squared goodness of fit statistic
chi2.chi_squared_exc = sum(((Ant_Sat.conc_exc - LSDAnt_chisquare.sat_exc).^2)./Ant_Sat.unc_exc.^2);
chi2.dof_exc = 4;

% Calculate the p-value
chi2.p_val_exc = 1 - chi2cdf(chi2.chi_squared_exc, chi2.dof_exc);

% Print the results 
fprintf('The chi-squared value is %.4f and the p-value is %.4f\n', chi2.chi_squared_exc, chi2.p_val_exc);

% This gives a p-value of 0.003 and a chi-squared value of 15 which is a
% better fit, but it is still not statistically significant. I don't know
% if there is a reason for these but I think if there was a geological
% explanation that could be a decent explanation for the fact that it is
% still a poor fit.










%% Plot of different production rate estimates 


% PR.est = [12.6, 12.0, 12.0, 12.1, 9.1, 12.6, 13.3, 13.7, 12.5, 15.84, 11.9, 13.9];
% PR.unc = [0.6, 1.3, 1.1, 0.8, 3.4, 1.6, 1.2, 1.7, 1.4, NaN, 1.4, 1.4];
% PR.pub = {'Bonneville Shoreline (UT, USA)', 'Scotland', 'New Zealand', 'Greenland', 'Arizona, USA', 'Death Valley (CA, USA)', 'White Mountains (CA, USA)', 'Chile', 'global', 'Theoretical', 'TU-Lab CRONUS-A','This Work'};
PR.est = [12.6, 12.0, 12.0, 12.1, 9.1, 12.6, 13.3, 13.7, 12.5, 15.84, 13.9];
PR.unc = [0.6, 1.3, 1.1, 0.8, 3.4, 1.6, 1.2, 1.7, 1.4, NaN, 1.4];
PR.pub = {'Bonneville Shoreline (UT, USA)', 'Scotland', 'New Zealand', 'Greenland', 'Arizona, USA', 'Death Valley (CA, USA)', 'White Mountains (CA, USA)', 'Chile', 'Global', 'Theoretical', 'This Work'};
PR.color = [0 0 0];


figure;

PR.x = 1:length(PR.pub);
errorbar(PR.x, PR.est, PR.unc, 'o', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineStyle', 'none', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'Color', 'k');

xticklabels(PR.pub);
xlim([0 12]);
ylim([8.5 16.5]);

xlabel('Calibration Dataset');
ylabel('Production Rate (atoms/g/yr)');

% Adjust figure aesthetics
set(gca, 'FontSize', 24, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5, 'XTick', 1:1:11); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure

%% Sensitivity of Equal Exposure Age to Production Rate

P14.refconc8ka = ((1 - exp(-8000.*const.lambda)).*LSD.P14SLHLmax)./const.lambda;
P14.exptime8ka = (-1./const.lambda).*log(1-(P14.refconc8ka.*const.lambda./P14.range));

P14.Pratio = P14.range/LSD.P14SLHLmax;

figure;
hold on
plot(P14.Pratio, P14.exptime8ka, 'k', 'LineWidth', 2);

% Define x and y values for peripheral targets
xvals = [0.9, 1, 1.1];
yvals = (-1./const.lambda) .* log(1 - (P14.refconc8ka .* const.lambda ./ (LSD.P14SLHLmax .* xvals)));

% Convert hex to RGB
c1 = sscanf('EFCEA6', '%2x%2x%2x', [1 3]) / 255;  % #EFCEA6
c2 = sscanf('B1A494', '%2x%2x%2x', [1 3]) / 255;  % #B1A494
c3 = sscanf('67553F', '%2x%2x%2x', [1 3]) / 255;  % #67553F

% Plot each point with its color
plot(xvals(1), yvals(1), '.', 'Color', c1, 'MarkerSize', 100)
plot(xvals(2), yvals(2), '.', 'Color', c2, 'MarkerSize', 100)
plot(xvals(3), yvals(3), '.', 'Color', c3, 'MarkerSize', 100)

    
xline(1, 'k--', 'LineWidth', 2);
yline(8000, 'k--', 'LineWidth',2);

text(0.9, 10600, 'A', 'FontSize', 54)
text(1.01, 8600, 'B', 'FontSize', 54)
text(1.11, 7350, 'C', 'FontSize', 54)


xlabel('Relative Production Rate (P/P_{true})')
ylabel('Exposure Age (years)')
yticklabels([4000:2000:22000])

% Adjust figure aesthetics
set(gca, 'FontSize', 54, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure

%xlim([0.85 1.15]);
%ylim([6000 10000]);

hold off


%% Sensitivity of UN-Equal Exposure Age to Production Rate
timerange = [200:200:26000];

for i=1:length(timerange)
    temp.refconc(i) = ((1 - exp(-1*timerange(i).*const.lambda)).*(LSD.P14SLHLmax+muon_PR))./const.lambda;
    temp.exptime(i,:) = (-1./const.lambda).*log(1-(temp.refconc(i).*const.lambda./P14.range));
    for j=1:10000
        if ~isreal(temp.exptime(i,j)) 
            temp.exptime(i,j) = NaN;
        end
    end
end

%%
figure;

P14.Pratio = P14.range/LSD.P14SLHLmax;

hold on
box on
cmap = (sky(length(timerange)));
colormap(sky(length(timerange)));
%colorbar();
% for i=1:length(timerange)
%     plot(P14.Pratio, temp.exptime(i,:), 'k', 'LineWidth', 2, 'Color', cmap(i, :));
% end

% Convert hex to RGB
c{1} = sscanf('DED4B9', '%2x%2x%2x', [1 3]) / 255;
c{2} = sscanf('BEA363', '%2x%2x%2x', [1 3]) / 255; 
c{3} = sscanf('A37C1C', '%2x%2x%2x', [1 3]) / 255;  

P14.ages = [6000, 12000, 18000];


yline(P14.ages(1), 'k:', 'LineWidth', 2)
yline(P14.ages(2), 'k:', 'LineWidth', 2)
yline(P14.ages(3), 'k:', 'LineWidth', 2)

for i=1:length(P14.ages)
    P14.ageconc(i) = ((1 - exp(-1*P14.ages(i).*const.lambda)).*LSD.P14SLHLmax)./const.lambda;
    P14.agecalc(i,:) = (-1./const.lambda).*log(1-(P14.ageconc(i).*const.lambda./P14.range));
    for j=1:10000
        if ~isreal(P14.agecalc(i,j)) 
            P14.agecalc(i,j) = NaN;
        end
    end
    plot(P14.Pratio, P14.agecalc(i,:), 'LineWidth', 4, 'Color', cell2mat(c(i)))
end


xline(1, 'k--', 'LineWidth', 2);

xlim([0.7 1.35]);
ylim([0 30000]);

xlabel('Relative Production Rate (P/P_{true})')
ylabel('Exposure Age (years)')
yticklabels([0:5000:30000])

% Adjust figure aesthetics
set(gca, 'FontSize', 54, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure


fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1200, 800];

%Turn on to save a high-resolution figure
%print(fig, 'prodratesensitivity.png', '-dpng', '-r600')

hold off

%% Paper figure related to the above 

P14.refconc20ka = ((1 - exp(-20000.*const.lambda)).*LSD.P14SLHLmax)./const.lambda;
P14.exptime20ka = (-1./const.lambda).*log(1-(P14.refconc20ka.*const.lambda./P14.range));

subplot(2,1,1)
plot(P14.Pratio, P14.exptime8ka);

subplot(2,1,2)
plot(P14.Pratio, P14.exptime20ka);




%%

elapsedTime = toc;
fprintf('Elapsed time to run this code: %.4f seconds\n', elapsedTime);
clear elapsedTime
close all
disp('Done. Have a great day :)')

%%%%%%%%%%%%%%%%%%%%%%%% End of Main Code %%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function matrix_data = createMatrix(extlab, conc)
    % createMatrix creates a matrix of lab data for plotting using the violin function.
    %
    % Inputs:
    %   - extlab: Cell array or array of labels identifying the lab
    %   - conc: Array of concentration values corresponding to each lab
    %   label
    %
    % Outputs:
    %   - matrix_data: A matrix where each column corresponds to a lab

    % Find unique labs
    unique_labs = unique(extlab);

    % Initialize a cell array to store the values for each lab
    split_data = cell(length(unique_labs), 1);

    % Split the data based on lab
    for i = 1:length(unique_labs)
        split_data{i} = conc(strcmp(extlab, unique_labs{i}))';
    end

    % Find the maximum length to pad the data to
    max_len = max(cellfun(@length, split_data));

    % Pad the vectors with NaNs to make them the same length
    padded_data = cellfun(@(x) [x, nan(1, max_len - length(x))], split_data, 'UniformOutput', false);

    % Combine the padded data into a matrix
    matrix_data = cell2mat(padded_data)';

end

function [lab, const] = structureLabData(extlab, conc, error, const)
    % structureLabData organizes data into a structured array for each lab.
    %
    % Inputs:
    %   - extlab: Cell array or array of labels identifying the lab for each conc.
    %   - conc: Array of concentration values.
    %   - error: Array of error values corresponding to each concentration.
    %   - const: Array of constants corresponding to the calculation
    %
    % Outputs:
    %   - lab: Structured array containing data for each lab.
    %   - const: Updated constant array with the unique lab information added 

    % Store a list of all the extraction labs in the constants
    const.labs = unique(extlab);

    % Initialize the structured array (1 for each lab with data)
    lab(1:length(const.labs)) = struct('name', '', 'conc', [], 'error', []);

    % Add data for each lab into a separate array in the structure
    for i = 1:length(const.labs)
        lname = const.labs{i}; % temp lab name variable
        idx = strcmp(extlab, lname); % get the indices for that specific lab
        lab(i).name = lname; % assign the extraction lab name to the array
        lab(i).conc = conc(idx); % add in the concentration measurements for that lab
        lab(i).error = error(idx); % add measurement errors for the concentration measurements
    end
end

function [St, lab] = StProdRate(site, name, lab, P14, const)
    % StProdRate determines the production rate of a sample scaled to SLHL
    % using the Stone, 2000 scaling framework
    %
    % Inputs:
    %   - site: relevant site data needed to perform the calculation
    %   (latitude, longitude, elevation, measured concentration).
    %   - name: Name of the site for ouput reference.
    %   - lab: The laboratory specific concentration information.
    %   - P14: The range of production rates to try calculated based on the
    %   literature (requires a P14.range vector of production rates over a
    %   range).
    %   - const: Constants assigned at the beginning of the code.
    %
    % Outputs:
    %   - St: A structured array containing all of the relevant information
    %   calculated during the production rate calibration in the function.
    %   - lab: save the laboratory specific production rate estimates in
    %   the lab structured array.
    %
    % The functions used within this code (stone2000 and ERA40) were
    % developed by Greg Balco, and unmodified for this purpose.
    %

    
    % Start with calculating the site specific (scaled up from SLHL input)
    % production rate for the site by calculating the scaling factor for
    % the site (term 2) and multiplying by the range of input production 
    % rates (term 1). Store the output in the St structured array.
    for i=1:length(P14.range)
        fieldname = ['P14z_', name]; % Make a more meaningful name
        St.(fieldname)(i) = P14.range(i) .* stone2000(site.lat, ERA40atm(site.lat, site.long, site.z), const.Fsp);
    end

    % Calculate the muon production rate to subtract from the total
    % reference production rate. Assumes average production rate of Balco
    % 2017 (3.07 atoms/g/yr) and Heisinger 2002 (3.78 atoms/g/yr) for total
    % muon SLHL surface production (uses 3.43 atoms/g/yr)
    St.muon = 3.43.*exp((1013.25-const.atmp)/const.L);

    
    % Calculate the saturation value due to spallation predicted for each production rate
    % includes a component of production due to muons here
    St.Nsat = (St.(fieldname)+St.muon)./(const.lambda);

    % Create a dictionary (key, value) where the key is the saturation
    % concentration at a given production rate and the value is the SLHL
    % scaled production rate that results in that saturation concentration.
    St.P14_map = dictionary(St.Nsat, P14.range);

    % Find the saturation value calculated from the above code that
    % minimizes the difference between the measured concentration and the
    % calculated concentration, then finds the corresponding value from the
    % dictionary.
    [~,avgidx] = min(abs(St.Nsat-site.avgdens));
    [~,maxidx] = min(abs(St.Nsat-site.maxlike));
    avgkey = St.Nsat(avgidx);
    maxkey = St.Nsat(maxidx);
    St.P14SLHLmax = St.P14_map(maxkey);
    St.P14SLHLavg = St.P14_map(avgkey);

    % Perform the same as above but for each lab
    for i=1:length(lab)
        [~,lab(i).avgidx] = min(abs(St.Nsat-lab(i).avg));
        [~,lab(i).maxidx] = min(abs(St.Nsat-lab(i).maxlike));
        lab(i).avgkey = St.Nsat(lab(i).avgidx);
        lab(i).maxkey = St.Nsat(lab(i).maxidx);
        lab(i).StP14SLHLmax = St.P14_map(lab(i).maxkey);
        lab(i).StP14SLHLavg = St.P14_map(lab(i).avgkey);
    end
end

function out = StProdRateUnc(site, P14, const)
    % StProdRateUnc determines the production rate of a sample scaled to SLHL
    % using the Stone, 2000 scaling framework without calculating the
    % laboratory production rates, it is a simplified version of the
    % function above intended to be used in Monte Carlo Analysis by varying
    % the inputs
    %
    % Inputs:
    %   - site: relevant site data needed to perform the calculation
    %   (latitude, longitude, elevation, measured concentration).
    %   - P14: The range of production rates to try calculated based on the
    %   literature (requires a P14.range vector of production rates over a
    %   range).
    %   - const: Constants assigned at the beginning of the code.
    %
    % Outputs:
    %   - out: The best fit production rate for each iteration of the monte-carlo.
    %
    % The functions used within this code (stone2000 and ERA40) were
    % developed by Greg Balco, and unmodified for this purpose.
    %

    
    % Start with calculating the site specific (scaled up from SLHL input)
    % production rate for the site by calculating the scaling factor for
    % the site (term 2) and multiplying by the range of input production 
    % rates (term 1). Store the output in the St structured array.
    for i=1:length(P14.range)
        PRrange(i) = P14.range(i) .* stone2000(site.lat, ERA40atm(site.lat, site.long, site.z), const.Fsp);
    end

    % Calculate the saturation value predicted for each production rate
    Nsat = (PRrange+site.muonPR)./(const.lambda);

    % Create a dictionary (key, value) where the key is the saturation
    % concentration at a given production rate and the value is the SLHL
    % scaled production rate that results in that saturation concentration.
    P14_map = dictionary(Nsat, P14.range);

    % Find the saturation value calculated from the above code that
    % minimizes the difference between the measured concentration and the
    % calculated concentration, then finds the corresponding value from the
    % dictionary.
    [~,idx] = min(abs(Nsat-site.conc));
    key = Nsat(idx);
    out = P14_map(key);
end

function out = LSDProdRateUnc(sampledata, flag, CA_conc, P14)
    % LSDProdRateUnc determines the production rate of a sample scaled to SLHL
    % using the Stone, 2000 scaling framework without calculating the
    % laboratory production rates, it is a simplified version of the
    % function above intended to be used in Monte Carlo Analysis by varying
    % the inputs
    %
    % Inputs:
    %   - site: relevant site data needed to perform the calculation
    %   (latitude, longitude, elevation, measured concentration).
    %   - P14: The range of production rates to try calculated based on the
    %   literature (requires a P14.range vector of production rates over a
    %   range).
    %   - const: Constants assigned at the beginning of the code.
    %
    % Outputs:
    %   - out: The best fit production rate for each iteration of the monte-carlo.
    %
    % The functions used within this code (stone2000 and ERA40) were
    % developed by Greg Balco, and unmodified for this purpose.
    %

    LSDnUncertainty = CD14C_CRONUScalib(sampledata, flag);
    

    temp.idx = find(LSDnUncertainty.tv<500000); % Find the indices of the time <50,000 years
    temp.sattv = LSDnUncertainty.tv(temp.idx); % Get the actual times for those indices

    % Next run through the output of the CD14C code to find the production rate
    % for each timestep corresponding to the trimmed time vector.
    for i=1:size(LSDnUncertainty.P14_CD, 1)
        temp.satP14(i,:) = LSDnUncertainty.P14_CD(i,temp.idx);
    end

    % Flip both vectors so that the accumulation starts 30 ka instead of in the
    % present.
    LSDnUncertainty.tv = flip(temp.sattv);
    LSDnUncertainty.P14v = fliplr(temp.satP14);


    % Eqution 4.1 of Dunai, 2010 assuming that the initial inventory C_inh is
    % zero calculates the saturation calculation for cosmogenic nuclides
    for b=1:size(LSDnUncertainty.P14v, 1)
        for a=1:length(LSDnUncertainty.tv)
            LSDnUncertainty.C14sum(b,a) = ((1 - exp(-1*LSDnUncertainty.tv(a).*sampledata.lambda)).*(LSDnUncertainty.P14v(b,a)+sampledata.muonPR))./sampledata.lambda;
        end
    end


    % Store the maximum value in a structured array for each access later
    max_fitconc = max((LSDnUncertainty.C14sum)');
    
    % If you change the P14.range somewhere else it will cause issues
    % because it is also calculated inside of the calibration function
    P14_map = dictionary(max_fitconc, P14);

    % Find the saturation value of the elevation scaled value that is closest 
    % to the CRONUS-A measurement from all of the compiled data and then 
    % determines what the SLHL scaled value is.
    [~,maxidxLSD] = min(abs(max_fitconc-CA_conc.conc));
    maxkeyLSD = max_fitconc(maxidxLSD);
    out = P14_map(maxkeyLSD);
end

function[h,L,MX,MED,bw]=violin(Y,varargin)

%defaults:
%_____________________
xL=[];
fc=[1 0.5 0];
lc='k';
alp=0.5;
mc='k';
medc='r';
b=[]; %bandwidth
plotlegend=1;
plotmean=1;
plotmedian=1;
x = [];
%_____________________

%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end

%get additional input parameters (varargin)
if isempty(find(strcmp(varargin,'xlabel')))==0
    xL = varargin{find(strcmp(varargin,'xlabel'))+1};
end
if isempty(find(strcmp(varargin,'facecolor')))==0
    fc = varargin{find(strcmp(varargin,'facecolor'))+1};
end
if isempty(find(strcmp(varargin,'edgecolor')))==0
    lc = varargin{find(strcmp(varargin,'edgecolor'))+1};
end
if isempty(find(strcmp(varargin,'facealpha')))==0
    alp = varargin{find(strcmp(varargin,'facealpha'))+1};
end
if isempty(find(strcmp(varargin,'mc')))==0
    if isempty(varargin{find(strcmp(varargin,'mc'))+1})==0
        mc = varargin{find(strcmp(varargin,'mc'))+1};
        plotmean = 1;
    else
        plotmean = 0;
    end
end
if isempty(find(strcmp(varargin,'medc')))==0
    if isempty(varargin{find(strcmp(varargin,'medc'))+1})==0
        medc = varargin{find(strcmp(varargin,'medc'))+1};
        plotmedian = 1;
    else
        plotmedian = 0;
    end
end
if isempty(find(strcmp(varargin,'bw')))==0
    b = varargin{find(strcmp(varargin,'bw'))+1}
    if length(b)==1
        disp(['same bandwidth bw = ',num2str(b),' used for all cols'])
        b=repmat(b,size(Y,2),1);
    elseif length(b)~=size(Y,2)
        warning('length(b)~=size(Y,2)')
        error('please provide only one bandwidth or an array of b with same length as columns in the data set')
    end
end
if isempty(find(strcmp(varargin,'plotlegend')))==0
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end
if isempty(find(strcmp(varargin,'x')))==0
    x = varargin{find(strcmp(varargin,'x'))+1};
end
%%
if size(fc,1)==1
    fc=repmat(fc,size(Y,2),1);
end

%% Calculate the kernel density
i=1;
for i=1:size(Y,2)
    
    if isempty(b)==0
        [f, u, bb]=ksdensity(Y{i},'bandwidth',b(i));
    elseif isempty(b)
        [f, u, bb]=ksdensity(Y{i});
    end
    
    f=f/max(f)*0.3; %normalize
    F(:,i)=f;
    U(:,i)=u;
    MED(:,i)=nanmedian(Y{i});
    MX(:,i)=nanmean(Y{i});
    bw(:,i)=bb;
    
end
%%
%-------------------------------------------------------------------------
% Put the figure automatically on a second monitor
% mp = get(0, 'MonitorPositions');
% set(gcf,'Color','w','Position',[mp(end,1)+50 mp(end,2)+50 800 600])
%-------------------------------------------------------------------------
%Check x-value options
if isempty(x)
    x = zeros(size(Y,2));
    setX = 0;
else
    setX = 1;
    if isempty(xL)==0
        disp('_________________________________________________________________')
        warning('Function is not designed for x-axis specification with string label')
        warning('when providing x, xlabel can be set later anyway')
        error('please provide either x or xlabel. not both.')
    end
end

%% Plot the violins
i=1;
for i=i:size(Y,2)
    if isempty(lc) == 1
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        end
    else
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        end
    end
    hold on
    if setX == 0
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)) ],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i)) ],[MED(:,i) MED(:,i)],medc,'LineWidth',2);
        end
    elseif setX == 1
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+x(i)-i],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+x(i)-i],[MED(:,i) MED(:,i)],medc,'LineWidth',2);
        end
    end
end

%% Add legend if requested
if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1
    
    if plotmean==1 & plotmedian==1
        L=legend([p(1) p(2)],'Mean','Median');
    elseif plotmean==0 & plotmedian==1
        L=legend([p(2)],'Median');
    elseif plotmean==1 & plotmedian==0
        L=legend([p(1)],'Mean');
    end
    
    set(L,'box','off','FontSize',32)
else
    L=[];
end

%% Set axis
if setX == 0
    axis([0.5 size(Y,2)+0.5, min(U(:)) max(U(:))]);
elseif setX == 1
    axis([min(x)-0.05*range(x) max(x)+0.05*range(x), min(U(:)) max(U(:))]);
end

%% Set x-labels
if isempty(xL)==0
    set(gca,'XTick',1:size(Y,2));  % Set x-axis ticks to match the number of violins
    set(gca,'XTickLabel',xL);      % Assign the x-axis labels directly from xL
end
% Hide only x-tick marks while keeping x-tick labels visible
ax = gca;
ax.XAxis.TickLength = [0 0];      % Remove x-axis tick marks
ax.YAxis.TickLength = [0.01, 0.025];    % Retain y-axis tick marks

% Adjust x-axis label angle and aesthetics
xtickangle(55);                   % Rotate x-axis labels by 55 degrees
set(gca, 'FontSize', 12, 'FontName', 'Helvetica');

% Keep y-axis ticks and box
box on;
%-------------------------------------------------------------------------
end %of function

function out = stone2000(lat,P, Fsp)

% Calculates the geographic scaling factor for cosmogenic-nuclide prodction as 
% a function of site latitude and atmospheric pressure, according to:
%
% Stone, J., 2000, Air Pressure and Cosmogenic Isotope Production. JGR 105:B10, 
% p. 23753. 
%
% Syntax: scalingfactor = stone2000(latitude,pressure,fsp)
%
% Units: 
% latitude in decimal degrees
% pressure in hPa
% fsp is the fraction (between 0 and 1) of production at sea level 
% and high latitude due to spallation (as opposed to muons). 
% This argument is optional and defaults to 0.978, which is the value 
% used by Stone (2000) for Be-10. The corresponding value for Al-26
% is 0.974. Note that using 0.844 for Be-10 and 0.826 for Al-26 will 
% closely reproduce the Lal, 1991 scaling factors as long as the standard
% atmosphere is used to convert sample elevation to atmospheric pressure.
% Also note that this function will yield the scaling factor for spallation
% only when fsp=1, and that for muons only when fsp=0.  
%
% Elevation can be converted to pressure with the functions
% stdatm.m (general use) and antatm.m (Antarctica). 
% 
% Vector argments are OK. All arguments must be the same size. 
%
% Written by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% First version, Feb. 2001
% checked March, 2006
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

% check for obvious errors

if ~isempty(find(abs(lat) > 90));
	error('Latitudes below 90, please');
end;

if length(lat) ~= length(P);
	error('Vectors the same size, please');
end;

% default Fsp

if nargin == 2;
	Fsp = 0.978;
end;

% Spallogenic production at index latitudes;

% enter constants from Table 1

a = [31.8518 34.3699 40.3153 42.0983 56.7733 69.0720 71.8733];
b = [250.3193 258.4759 308.9894 512.6857 649.1343 832.4566 863.1927];
c = [-0.083393 -0.089807 -0.106248 -0.120551 -0.160859 -0.199252 -0.207069];
d = [7.4260e-5 7.9457e-5 9.4508e-5 1.1752e-4 1.5463e-4 1.9391e-4 2.0127e-4];
e = [-2.2397e-8 -2.3697e-8 -2.8234e-8 -3.8809e-8 -5.0330e-8 -6.3653e-8 -6.6043e-8];

ilats = [0 10 20 30 40 50 60];

% calculate index latitudes at given P's

lat0 = a(1) + (b(1) .* exp(P./(-150))) + (c(1).*P) + (d(1).*(P.^2)) + (e(1).*(P.^3));
lat10 = a(2) + (b(2) .* exp(P./(-150))) + (c(2).*P) + (d(2).*(P.^2)) + (e(2).*(P.^3));
lat20 = a(3) + (b(3) .* exp(P./(-150))) + (c(3).*P) + (d(3).*(P.^2)) + (e(3).*(P.^3));
lat30 = a(4) + (b(4) .* exp(P./(-150))) + (c(4).*P) + (d(4).*(P.^2)) + (e(4).*(P.^3));
lat40 = a(5) + (b(5) .* exp(P./(-150))) + (c(5).*P) + (d(5).*(P.^2)) + (e(5).*(P.^3));
lat50 = a(6) + (b(6) .* exp(P./(-150))) + (c(6).*P) + (d(6).*(P.^2)) + (e(6).*(P.^3));
lat60 = a(7) + (b(7) .* exp(P./(-150))) + (c(7).*P) + (d(7).*(P.^2)) + (e(7).*(P.^3));

% initialize output

correction = zeros(size(P));

% northernize southern-hemisphere inputs

lat = abs(lat);

% set high lats to 60;

lat(find(lat > 60)) = (zeros(size(find(lat > 60))) + 60);

% loop 

b =1;

while b <= length(lat);		

	%interpolate for actual elevation:

	S(b) = interp1(ilats,[lat0(b) lat10(b) lat20(b) lat30(b) lat40(b) lat50(b) lat60(b)], lat(b));
	
	% continue loop	

	b = b+1;
	
end;

% Production by muons

% constants

mk = [0.587 0.600 0.678 0.833 0.933 1.000 1.000];

% index latitudes at given P's

ml0 = mk(1) .* exp((1013.25 - P)./242);
ml10 = mk(2) .* exp((1013.25 - P)./242);
ml20 = mk(3) .* exp((1013.25 - P)./242);
ml30 = mk(4) .* exp((1013.25 - P)./242);
ml40 = mk(5) .* exp((1013.25 - P)./242);
ml50 = mk(6) .* exp((1013.25 - P)./242);
ml60 = mk(7) .* exp((1013.25 - P)./242);

% loop 

b =1;

while b <= length(lat);		

	%interpolate for actual elevation:

	M(b) = interp1(ilats,[ml0(b) ml10(b) ml20(b) ml30(b) ml40(b) ml50(b) ml60(b)], lat(b));
	
	% continue loop	

	b = b+1;
	
end;

% Combine spallogenic and muogenic production; return

Fm = 1 - Fsp;

out_1 = ((S .* Fsp) + (M .* Fm));

% make vectors horizontal

if size(out_1,1) > size(out_1,2);
	out = out_1';
else;
	out = out_1;
end;
end;
