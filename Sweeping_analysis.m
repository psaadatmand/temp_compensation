
% To calculate the scaling factor for different temperature range change
% the "n" in the following line
% bin_edges = 10 : n : 45

% Load the data from Excel file
% Replace 'DATA_7%.xlsx' with your actual file name if needed.
data = readtable('DATA_8%.xlsx');


% Extract the temperature and tau value columns
temperature = data.Temperature;
tau = data.Tau_value;

% Define the temperature bin edges (from 10°C to 45°C, with 1 degree intervals)
bin_edges = 10:2:45;

% Bin the temperature data into the defined bins
[~, ~, bin_indices] = histcounts(temperature, bin_edges);

% Initialize an array to store the mean tau values for each bin
mean_tau_per_bin = NaN(length(bin_edges)-1, 1);
std_tau_per_bin = NaN(length(bin_edges)-1, 1);  % To store the standard deviation


% Calculate the mean tau value for each bin
for i = 1:length(mean_tau_per_bin)
    tau_in_bin = tau(bin_indices == i);  % Tau values within the current bin
    if ~isempty(tau_in_bin)
        mean_tau_per_bin(i) = mean(tau_in_bin);  % Calculate mean for the current bin
        std_tau_per_bin(i) = std(tau_in_bin);    % Standard deviation for the current bin

    end
end

% Create temperature ranges as strings
temperature_ranges = cell(length(bin_edges)-1, 1);
for i = 1:length(temperature_ranges)
    temperature_ranges{i} = sprintf('%d-%d', bin_edges(i), bin_edges(i+1));
end

% Display the results
disp('Temperature Range (°C)   Mean Tau Value');
for i = 1:length(mean_tau_per_bin)
    fprintf('%2d - %2d°C:   %f\n', bin_edges(i), bin_edges(i+1), mean_tau_per_bin(i));
end

% Calculate the midpoint of each temperature bin for plotting
temp_bin_midpoints = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;


% Find the mean tau for the 22°C to 23°C temperature range
index_22 = find(temp_bin_midpoints >= 22 & temp_bin_midpoints < 25);
mean_tau_22 = mean_tau_per_bin(index_22);



% Calculate the fraction mean_tau / mean_tau(temperature 22-23)
scalingfactor_tau = mean_tau_22 ./mean_tau_per_bin;
error_tau = (mean_tau_per_bin - mean_tau_22)/mean_tau_22;



% Plot the fraction vs temperature
figure;
plot(temp_bin_midpoints, scalingfactor_tau, '--o', 'LineWidth', 1);
xlabel('Temperature (°C)');
ylabel('Scaling factor of Mean Tau');
title('Scaling factor of Mean Tau vs Temperature O2 = 7% (Normalized to 22-23°C)');
grid on;

% Set the grid and axis ticks
xticks(10:2:45);  % Set x-axis ticks from 10°C to 45°C with 1-degree intervals
ax = gca;
ax.XGrid = 'on';   % Turn on vertical grid lines for x-axis
ax.YGrid = 'on';   % Turn on horizontal grid lines for y-axis

hold off;


% Plot the fraction vs temperature
figure;
plot(temp_bin_midpoints, error_tau, 'r--o', 'LineWidth', 1);
xlabel('Temperature (°C)');
ylabel('Fraction of Mean Tau');
title('Error of Mean Tau vs Temperature O2 = 7% (Normalized to 22-23°C)');
grid on;

% Set the grid and axis ticks
xticks(10:2:45);  % Set x-axis ticks from 10°C to 45°C with 1-degree intervals
ax = gca;
ax.XGrid = 'on';   % Turn on vertical grid lines for x-axis
ax.YGrid = 'on';   % Turn on horizontal grid lines for y-axis

hold off;

% Create a table with Temperature Midpoints, Mean Tau, and Fraction Tau
scaling_factor_table= table(temperature_ranges , mean_tau_per_bin, scalingfactor_tau, ...
    'VariableNames', {'Temperature range', 'Mean_Tau', 'scaling factor'});

Error_table = table(temperature_ranges , mean_tau_per_bin, error_tau, ...
    'VariableNames', {'Temperature range', 'Mean_Tau', 'Error_Tau'});

% Display the table
disp(scaling_factor_table);
disp(Error_table);

% Plot the mean tau values against the temperature bin midpoints
figure;
hold on;
plot(temp_bin_midpoints, mean_tau_per_bin, 'g-o', 'LineWidth', 2);
plot(temperature, tau,'.r','LineWidth',2);
xlabel('Temperature (°C)');
ylabel('Mean Tau Value');
title('Mean Tau Value vs Temperature in  O2 = 9%');
hold off;

xticks(10:2:45);
grid on;
ax = gca;
ax.XGrid = 'on';   % Turn on vertical grid lines for x-axis
ax.YGrid = 'on'; 

% Plot the mean tau values against the temperature bin midpoints
figure;
hold on;
plot(temp_bin_midpoints, mean_tau_per_bin, 'r-o', 'LineWidth', 2);
plot(temp_bin_midpoints, calibrated_mean_tau, 'g-o', 'LineWidth', 2);
xlabel('Temperature (°C)');
ylabel('Mean Tau Value');
title('Mean Tau Value vs Temperature in  O2 = 9%');
hold off;

xticks(10:2:45);
grid on;
ax = gca;
ax.XGrid = 'on';   % Turn on vertical grid lines for x-axis
ax.YGrid = 'on'; 


%plot Error bar
figure;
hold on;
errorbar(temp_bin_midpoints, mean_tau_per_bin, std_tau_per_bin, 'g-o','LineWidth', 2);
xlabel('Temperature (°C)');
ylabel('Mean Tau Value');
title('Mean Tau Value vs Temperature in  O2 = 7%');
hold off;


% Set the grid and axis ticks
xticks(10:2:45);
grid on;
ax = gca;
ax.XGrid = 'on';   % Turn on vertical grid lines for x-axis
ax.YGrid = 'on';   % Turn on horizontal grid lines for y-axis








