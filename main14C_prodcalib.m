% This code was originally written by Keir Nichols as a PhD student at
% Tulane University to calculate the spallation production rate of in-situ 14C in
% quartz using the CRONUS-A measurement as a calibration measurement. It was
% subsequently updated by Jason Drebber as a PhD Student at the Colorado
% School of Mines to include a larger number of CRONUS-A measurements from
% different labs and expanded to incorporate the LSDn scaling framework and
% more statistical tests


% clear figures, output in command window and workspace variables
clf, close, clc, clear

%% Define Constants and Load Data
const.minlambda = log(2)/5670; % 14C decay constants (calculated using a half-life of
const.lambda = log(2)/5700;    % of 5700 +/- 30 years following Hippe and Lifton, 2014;
const.maxlambda = log(2)/5730; % This is the IAEA accepted value as of June 2024)

const.Fsp = 1; % Used in stone scaling to account for muogenic production (Balco, 2008)

CA.lat = -77.883; % Latitude of CRONUS-A sample (Jull et al., 2015)
CA.long = 160.9431; % Longitude of CRONUS-A sample (Jull et al., 2015)
CA.z = 1612; % Elevation of CRONUS-A site in meters (Jull et al., 2015)
CA.z_error = 10; % Error in elevation measurements from a GPS collected point (Lecavlier, 2022)

% Load data
load CRONUSA.mat % CRONUS-A measurements
load extlab.mat % Extraction Lab data
CA.conc = CRONUSA(:,4); % CA is a structure of CRONUS-A data; conc in the concentration measurements
CA.error = CRONUSA(:,5); % Published error on the concentration measurements


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

%Fix this so its not hard coded   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const.labs = {"Purdue", "P-CEGS", 'TU-CEGS', 'LDEO', 'ANST0-UOW', 'ETH-Zurich', 'Intercomparison', 'University of Cologne'};

% Reorder the data to correspond to the lab name order
lab = lab([6, 5, 7, 4, 1, 2, 3, 8]);

% Store the lab data without the University of Cologne
labtrim = lab(1:7);

% Remove the University of Cologne data from the compiled data because
% there is low certainty in the values
[valid_idx] = find(extlab ~= "University of Cologne");
CA.conc = CA.conc(valid_idx);
CA.error = CA.error(valid_idx);
extlab = extlab(valid_idx);

clear valid_idx


%% Plot the main compilation figure
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

CA.tot_density = temp.tot_density;

% Add axis labels with refined font sizes
xlabel('^{14}C concentration (10^{5} atoms g^{-1})', 'FontSize', 16, 'FontWeight', 'Bold', 'FontName', 'Helvetica');
ylabel('Probability Density', 'FontSize', 16, 'FontWeight', 'Bold', 'FontName', 'Helvetica');

ax = gca;  % Get current axes
ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
ax.XAxis.Exponent = 0;
ax.XAxis.TickLabelFormat = '%.1f'; % Keep scientific notation
ax.YTickLabel = [];
ax.YTick = [];

% Define legend categories and colors
temp.unique_categories = {'Purdue', 'P-CEGS', 'TU-CEGS', 'LDEO', ...
                          'ANSTO-UOW', 'ETH Zurich', 'Intercomparison'};
legend_colors = [0.121, 0.467, 0.706;   % Purdue
                 0.172, 0.627, 0.745;   % P-CEGS
                 0.850, 0.372, 0.007;   % TU-CEGS
                 0.702, 0.871, 0.412;   % LDEO
                 0.525, 0.396, 0.750;   % ANSTO-UOW
                 0.984, 0.502, 0.447;   % ETH Zurich
                 0.596, 0.596, 0.596];  % Intercomparison

% Initialize dummy patch handles for legend
temp.legend_handles = gobjects(1, length(temp.unique_categories));

% Create transparent square patches for legend
for i = 1:length(temp.unique_categories)
    temp.legend_handles(i) = patch(NaN, NaN, legend_colors(i, :), ...
        'FaceAlpha', 0.6, 'EdgeColor', 'none'); % Transparent square
end

% Generate the legend with transparent patches
leg = legend(temp.legend_handles, temp.unique_categories, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'FontSize', 16, 'Box', 'off', 'TextColor', 'k', 'FontName', 'Helvetica');

% Force legend markers to be square
leg.ItemTokenSize = [12, 12]; % Set size for legend markers (width, height)

% Set consistent axis limits
xlim([4.5 8.5]); % Use a fixed range for x-axis
ylim([0,max(temp.tot_scale_density)]); % Add some buffer for the y-axis

box on; % Adds a box around the current axes

% Adjust figure aesthetics
set(gca, 'FontSize', 16, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure

fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1200, 800]; % Example: 1200x800 pixels

%Turn on to save a high-resolution figure
%print(fig, 'CombinedDensity.png', '-dpng', '-r600')

hold off;
clear lab_name temp i ax legend_colors leg fig

%% Calculate the density of the adjusted data
% This section will calculate the adjusted lab kernel density which means
% the total density without the Intercomparison material. This is mean to
% provide a baseline comparison for the new data  to compare to the data
% which includes the intercomparison material (calculated above in the
% figure producing section.

% Remove Intercomparison data from the dataset
idxlab = ~strcmp(extlab, "Intercomparison");
adjlab = extlab(idxlab);
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


%% Plot the individual lab data each on their own plot to compare
% This figure will plot the cumulative kernel density for each lab on its
% own subplot which allows for easier identification of trends within each
% lab. It also shows a violin plot (similar to a boxplot) of the data which
% shows how the mean of each lab compares to the mean of each other lab.

% Make a new figure
figure;

% Subplot positions defined manually for figure structure
subplot_positions = [
    0.09, 0.68, 0.25, 0.24;  % Top-left
    0.34, 0.68, 0.25, 0.24;  % Top-center
    0.59, 0.68, 0.25, 0.24;  % Top-right
    0.09, 0.44, 0.25, 0.24;  % Middle-left
    0.34, 0.44, 0.25, 0.24;  % Middle-center
    0.59, 0.44, 0.25, 0.24;  % Middle-right
    0.09, 0.20, 0.25, 0.24;  % Bottom-left
    0.34, 0.20, 0.25, 0.24;  % Bottom-center
    0.59, 0.20, 0.25, 0.24;  % Bottom-right
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
    text(0.05, 0.9, [temp.subcaptions{i} lab(i).name], 'Units', 'normalized', "FontSize", 16, 'FontWeight', 'normal', 'FontName', 'Helvetica');
    set(gca, 'FontSize', 16, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis

    % Add axis labels and adjust positioning conditionally
    if i == 1 || i == 4
        ylabel("Probability Density");
        set(gca, 'XTickLabel', [])
    elseif i == 8
        xlabel('^{14}C Concentration (10^{5} atoms/g)');
        set(gca, 'YTickLabel', [])
    elseif i == 7
        xlabel('^{14}C Concentration (10^{5} atoms/g)'); 
        ylabel("Probability Density");
    elseif i == 3
        xlabel('^{14}C Concentration (10^{5} atoms/g)'); 
        set(gca, 'XAxisLocation', 'top');
    end
    
    % Remove axis labels for specific subplots
    if ismember(i, [2, 5, 6])
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
    end
end

% Add violin plots in the bottom right with adjusted position
axes_handles(end) = subplot('Position', subplot_positions(9, :));
violin(CA.matrix_data/1e5, 'xlabel', cellstr(const.labs), 'facecolor', [1 1 1], 'medc', []);
ylabel("^{14}C Concentration (10^{5} atoms/g)");
text(0.05, 0.9, 'I. Lab Comparison', 'Units', 'normalized', "FontSize", 16, 'FontWeight', 'normal', 'FontName', 'Helvetica');
yticks(5:1:10);
set(gca, 'YAxisLocation', 'right');  % Move x-axis labels to the top

% Link axes for easy comparison
linkaxes(axes_handles(1:end-1), 'xy');

% Set some specific features for the figure to make it look nicer
set(gca, 'FontSize', 16, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure

% Set figure size
fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1200, 800]; % Example: 1200x800 pixels

%Turn on to save a high-resolution figure
%print(fig, 'labspecific.png', '-dpng', '-r600')

hold off;
clear i idx lab_name temp ax fig subplot_positions ans axes_handles;


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


%% Establish Production Rate Vector
% Establish a range of SLHL production rate values to try based on
% an estimated range from the literature.

% To extend the range of values or decrease the spacing (default is 0.001
% atoms/g/yr) modify the min, max and interval variables below.

P14.min = 10; %minimum estimate
P14.max = 20; %maximum estimate
P14.intervals = (P14.max-P14.min)*1000; %number of interval values to check
P14.range = linspace(P14.min, P14.max, P14.intervals); %vector of production rate estimates
 

%% Stone Scaling Production Rate Estimate
% This code block runs the StPRodRate function which performs a production
% rate calibration of an input site (in this case the CRONUS-A site) and
% calibrates the production rate for each lab where there is data using the
% Stone, 2000 scaling scheme to scale the estimate to SLHL and the ERA40
% atmospheric model. 


[St, lab] = StProdRate(CA, 'CRONUSA', lab, P14, const);

%Isochrons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t=2000:2000:20000;  %plots isochrons for 2000 kyr exposure durations up to 20 kyr
% for i=1:length(t)
%     N=(P14z/lambda)*(1-exp(-t*lambda));
% end

%% Monte-Carlo estimate of CRONUS-A distribution
% This section performs a monte-carlo analysis of the 

for i=1:10000
    while true
        % Draw a random value with uniform probability over the full range of
        % the CRONUS-A data
        temp.x_rand(i) = min(conc_range.full_comp) + (max(conc_range.full_comp)-min(conc_range.full_comp))*rand;

        % Calculate the PDF probability for that value using a linear
        % interpolation
        temp.y_rand(i) = interp1(conc_range.full_comp, CA.tot_density, temp.x_rand(i), 'linear', 0);

        % Use a conditional to accept or not accept the value based on the pdf
        % distribution
        if rand < temp.y_rand(i)
            CA_mc(i) = temp.x_rand(i);
            break;
        end
    end
end

ksdensity(CA_mc)
%%

 for i=1:length(P14.range)
        test(i) = P14.range(i) .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, CA.z), const.Fsp);
 end

 test_sat = test./(const.lambda);

 test_map = dictionary(test_sat, P14.range);

 for i=1:length(CA_mc)
    [~,idx] = min(abs(test_sat-CA_mc(i)));
    test_key = St.Nsat(idx);
    test_PR(i) = St.P14_map(test_key);
 end

 %%
 test_PR = sort(test_PR, 'ascend');
 test_dens = kde(test_PR, "NumPoints", 5000);
 [~,I] = max(test_dens)
 test_PR(I)


%% LSDn Scaling Section
% This section calculates the time dependentproduction rate using the LSDn 
% Scaling (Lifton et al., 2014). The code here is from the CD14C code published
% in Koester and Lifton, 2023 which calculates the compositionally
% dependent production rate using a P14 of 13.5 for quartz from
% spallation. This version is modified to take a range of possible
% production values (corresponding to the same range used in the St scaling
% method) to output a range of time dependent production rates that can be
% used to determine the best fit production rate.
    

LSD = CD14C_CRONUScalib('LSDn_inputs.txt', 1);
% The output here "LSD" is a structure that stores the time and
% compostitionally dependent production rate at specific time intervals
% given by the vector "tv" defining the time period of each production rate
% estimate.


%% Calculate the saturation concentration
% First find the indices of the time vector less than 50,000 years. We expect any
% continuously exposed material that is not undergoing erosion to be satured
% with respect to in-situ 14C after around 30,000 years, so 50,000 years is
% a conservative estimate. This step reduces processing time.

temp.idx = find(LSD.tv<500000); % Find the indices of the time <50,000 years
temp.sattv = LSD.tv(temp.idx); % Get the actual times for those indices

% Next run through the output of the CD14C code to find the production rate
% for each timestep corresponding to the trimmed time vector.
for i=1:size(LSD.P14_CD, 1)
    temp.satP14(i,:) = LSD.P14_CD(i,temp.idx);
end

% Flip both vectors so that the accumulation starts 30 ka instead of in the
% present.
LSD.tv = flip(temp.sattv);
LSD.P14v = fliplr(temp.satP14);


% Eqution 4.1 of Dunai, 2010 assuming that the initial inventory C_inh is
% zero calculates the saturation calculation for cosmogenic nuclides
for b=1:size(LSD.P14v, 1)
    for a=1:length(LSD.tv)
        LSD.C14sum(b,a) = ((1 - exp(-1*LSD.tv(a).*const.lambda)).*LSD.P14v(b,a))./const.lambda;
    end
end


% Store the maximum value in a structured array for each access later
LSD.max = max((LSD.C14sum)');

clear temp;

%% Plot the time-dependent accumulation
% This plot shows the accumulation of cosmogenic nuclides up to the
% saturation concentration using the equation 4.1 from Dunai, 2010. The
% plot is meant to demonstrate the relationship.

% Plot the log log relationship because the features are on the log scale
figure;
loglog(LSD.tv, LSD.C14sum(4705,:), 'LineWidth', 2,...
       'DisplayName', '^{14}C Concentration Curve');
grid on; % Add gridlines

% Customize the grid lines to make them less bright
ax = gca; % Get the current axes
ax.GridAlpha = 0.3; % Set grid transparency (lower value for less brightness)
ax.GridColor = [0.5, 0.5, 0.5]; % Set grid color (lighter gray)

% Customize axes labels and make them larger
xlabel('Exposure Duration (years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('In-situ ^{14}C Concentration (atoms/g)', 'FontSize', 14, 'FontWeight', 'bold');

% Add a horizontal line at the saturation concentration
saturation_concentration = round(LSD.max(4705)); % Get the saturation concentration as an integer
hline = yline(saturation_concentration, '--r', 'LineWidth', 2, ...
    'DisplayName', 'Saturation Concentration'); % Dashed red line for visibility

% Add text with customized font size and move it to align with the legend
text(5900, 35000, sprintf('Saturation: %d', saturation_concentration), ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', 'black'); % Original position

% Add a legend with both the curve and horizontal line
legend('show', 'Location', 'best', 'FontSize', 12); % Align legend near the text

% Customize tick marks
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Larger axis ticks and thicker axis lines


%% Estimate the production rate
% This section will use the saturation concentrations that were calculated
% above to establish a dictionary of saturation values and production rates
% that determine those values, then it will calculate the best fit of the
% production rate based on the saturation values.

% Create a dictionary of the saturation values and the reference production
% rates that result in those saturation values
LSD.P14_map = dictionary(LSD.max, P14.range);

% Find the saturation value of the elevation scaled value that is closest 
% to the CRONUS-A measurement from all of the compiled data and then 
% determines what the SLHL scaled value is.
[~,avgidxLSD] = min(abs(LSD.max-CA.avgdens));
[~,maxidxLSD] = min(abs(LSD.max-CA.maxlike));
avgkeyLSD = LSD.max(avgidxLSD);
maxkeyLSD = LSD.max(maxidxLSD);
LSD.P14SLHLmax = LSD.P14_map(maxkeyLSD);
LSD.P14SLHLavg = LSD.P14_map(avgkeyLSD);

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

%% Calculate the elevation scaled saturation curve for different estimates
% This section calculates the saturation curve using Stone scaling to
% estimate the saturation value dependent on elevation. 

z_range=0:100:3000; %elevation range to calculate the values

% For loop that calculates the scaled production rate at elevation based 
% on estimated SLHL production rates calculated using different statistical
% likelihoods.
for i=1:length(z_range)
    St.P14avg(i,1) = St.P14SLHLavg .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp);
    St.P14max(i,1) = St.P14SLHLmax .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp);
    for j=1:length(lab)
        lab(j).P14avg(i,1) = lab(j).StP14SLHLavg .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp)
        lab(j).P14max(i,1) = lab(j).StP14SLHLmax .* stone2000(CA.lat, ERA40atm(CA.lat, CA.long, z_range(i)), const.Fsp)
    end
end

% Calculate the saturation curve concentration for different elevations
% scaled based on the production rate at those elevations.
St.Nsatavg = St.P14avg/(const.lambda)
St.Nsatmax = St.P14max/(const.lambda)

for i=1:length(lab)
    lab(i).Nsatavg = lab(i).P14avg/(const.lambda)
    lab(i).Nsatmax = lab(i).P14max/(const.lambda)
end

%% Calculate the saturation curve for the LSD scaling framework
LSDsatcurvein = {};

for i=1:length(z_range)
    LSDsatcurvein{i,1} = 'CRONUS-A';
    LSDsatcurvein{i,2} = CA.lat;
    LSDsatcurvein{i,3} = CA.long;
    LSDsatcurvein{i,4} = z_range(i);
    LSDsatcurvein{i,5} = 100;
    for ii=1:10
        LSDsatcurvein{i,ii+5} = 0;
    end
end

satcurve_test = CD14C_CRONUScalib(LSDsatcurvein, 2, LSD.P14SLHLmax);


%% Calculate the saturation values for each elevation

idx = find(satcurve_test.tv<500000);
sattv = satcurve_test.tv(idx); %Find the values

for i=1:size(satcurve_test.P14_CD, 1)
    satcurveP14(i,:) = satcurve_test.P14_CD(i,idx); %Find the production rates for the time
end

%Flip both vectors so that the accumulation starts 30 ka instead of in the
%present.
satcurve_test.tvt = flip(sattv);
satcurve_test.P14v = fliplr(satcurveP14);

for b=1:size(satcurve_test.P14v, 1)
    for a=1:length(satcurve_test.tvt)
        LSDsatcurve.C14sum(b,a) = ((1 - exp(-1*satcurve_test.tvt(a).*const.lambda)).*satcurve_test.P14v(b,a))./const.lambda;
    end
end

LSD.maxcurve = max((LSDsatcurve.C14sum)');


%% Plot the saturation curves and all Antarctic 14C data from ICE-D

Ant_Sat.conc = [183030; 968970; 160050; 974370; 1177930; 1038010];
Ant_Sat.lat = [-70.86; -77.75; -70.82; -77.75; -73.44; -73.39];
Ant_Sat.long = [68.13; 160.8; 68.17; 160.8; 61.9; 61.72];
Ant_Sat.elev = [225; 2160; 100; 2020; 2538; 2137];
Ant_Sat.unc = [8420; 15770; 12860; 19180; 19490; 20640];
Ant_Sat.names = {"98-PCM-010-SRDK"; "WBC-UVP"; "98-PCM-002-BVLK"; "WBC-2020"; "98-PCM-105-MNZ"; "98-PCM-067-MNZ"};



load all_Antarctic_14C.txt %text file with all in situ 14C measurements in ICE-D in Antarctica (as of 3 May 2024)
data1=all_Antarctic_14C;
ant.z=data1(:,1); ant.conc=data1(:,2); err=data1(:,3); sixpercenterr=data1(:,4);  %column 3 = analytical uncertainties, 4 = 6%



figure;
hold on
St.curvemax = plot(St.Nsatmax/1e5,z_range);
LSD.curvemax = plot(LSD.maxcurve/1e5, z_range);
%St.punc = plot(Nsat2,z);
%St.nunc = plot(Nsat3,z);
%set(St.curveavg, 'LineWidth',3);
%set(St.curveavg, 'color',[0 0 0]);
set(St.curvemax, 'LineWidth',3, 'color',[0.937, 0.502, 0.502]);
%set(St.curvemax, 'color',[0.937, 0.502, 0.502]);
set(LSD.curvemax, 'LineWidth',4);
set(LSD.curvemax, 'color',[0.529, 0.808, 0.980]);
% set(h2, 'LineWidth',2)
% set(h2, 'color',[0.9290 0.6940 0.1250])
% set(h3, 'LineWidth',2)
% set(h3, 'color',[0.9290 0.6940 0.1250])

%plot(N,z,':k') %plots isochrons

%Plot the CRONUS-A value
plot(CA.maxlike/1e5,CA.z,'o', MarkerSize= 12, markerfacecolor =[0.882, 0.745, 0.416], MarkerEdgeColor='k')

%Plot the other Antarctic Saturated Surfaces
testsamples = plot(Ant_Sat.conc/1e5, Ant_Sat.elev, 'o', MarkerSize= 10, markerfacecolor = [0.251, 0.690, 0.651], MarkerEdgeColor='k');

ant.curve = plot(ant.conc/1e5,ant.z,'o');    % plots all in situ 14C measurements with no error bars
%h5 = errorbar_x(samples,elev1,err,'square')     % above but with analytical uncertainties
%h5 = errorbar_x(samples,elev1,sixpercenterr,'square')  %above but 6 % uncertainties
set(ant.curve, 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0], 'Marker','.');


%h4 = errorbar_x(samples3,elev3,cronus6percent,'.')
%set(h4, 'MarkerSize', 20, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0 0 0])

xlabel('^{14}C concentration (10^{5} atoms g^{-1})')
ylabel('Elevation (m asl)')
xlim([0 15])
ylim([0 3000])
legend('St Saturation Curve', 'LSDn Saturation Curve', 'CRONUS-A Maximum Likelihood', 'Other Saturated Surfaces', 'Antarctic Data',  'Box', 'off')
legend('Position', [0.65, 0.3, 0.1, 0.1]);
set(gca, 'XTick', [0:2.5:15])

hold off;
box on;

% Adjust figure aesthetics
set(gca, 'FontSize', 16, 'FontWeight', "Normal", 'FontName', 'Helvetica', 'LineWidth', 1.5); % Adjust tick labels and axis
set(gcf, 'Color', 'w'); % White background for the figure

fig = gcf; % Get current figure handle
fig.Position = [100, 100, 1200, 800]; % Example: 1200x800 pixels

%Turn on to save a high-resolution figure
%print(fig, 'saturationcurve.png', '-dpng', '-r600')


%% Calculate the saturation curve for the LSD scaling framework
Sat_val = {};

for i=1:length(Ant_Sat.elev)
    Sat_val{i,1} = Ant_Sat.names(i);
    Sat_val{i,2} = Ant_Sat.lat(i);
    Sat_val{i,3} = Ant_Sat.long(i);
    Sat_val{i,4} = Ant_Sat.elev(i);
    Sat_val{i,5} = 100;
    for ii=1:10
        Sat_val{i,ii+5} = 0;
    end
end


other_Ant = CD14C_CRONUScalib(Sat_val, 2, LSD.P14SLHLmax);


%% Calculate the saturation values for each test site in Antarctica for the CHi-squared test

idx = find(other_Ant.tv<500000);
sattv = other_Ant.tv(idx); %Find the values

for i=1:size(other_Ant.P14_CD, 1)
    Ant_chisquare(i,:) = other_Ant.P14_CD(i,idx); %Find the production rates for the time
end

%Flip both vectors so that the accumulation starts 30 ka instead of in the
%present.
other_Ant.tvt = flip(sattv);
other_Ant.P14v = fliplr(Ant_chisquare);

for b=1:size(other_Ant.P14v, 1)
    for a=1:length(other_Ant.tvt)
        LSDAnt_chisquare.C14sum(b,a) = ((1 - exp(-1*other_Ant.tvt(a).*const.lambda)).*other_Ant.P14v(b,a))./const.lambda;
    end
end

LSDAnt_chisquare.sat = max((LSDAnt_chisquare.C14sum)')';

%% Calculate chi squared

chi_squared = sum(((Ant_Sat.conc - LSDAnt_chisquare.sat).^2)./Ant_Sat.unc.^2)
dof = 6;

p_val = 1 - chi2cdf(chi_squared, dof)

%Wlep this gives a p-value of 0 and a chi-squared value of 168 which is a
%bad fit, but there are two points that plot clearly below saturation,
%which tells me that maybe they are causing this because the other 4
%samples are at saturation. 


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


    %We may need to take into account an uncertainty on the elevation
    %measurement, too - CRONUS-A was collected by Greg in 2004
    %the 0.8 value at the end is the fraction of production at SLHL due to spallation (as
    %opposed to muons) - I've used 0.8, but we might need to vary this. E.g.
    %Greg mentions "~15 %" and "20 %" in his 2017 paper
    %New thought, we should probably provide a range of estimates based on the
    %ful range of values reported in the literature from Dyonisius to Balco.

    % Calculate the saturation value predicted for each production rate
    St.Nsat = St.(fieldname)./(const.lambda);


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
    
    set(L,'box','off','FontSize',14)
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
