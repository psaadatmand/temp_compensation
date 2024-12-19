% Load your data into MATLAB
% Assuming the data is stored in a table called 'data'
% with columns 'Tau_value' and 'Temperature'

data_O2_6 = readtable('‌DATA_6%.xlsx');
data_O2_7 = readtable('DATA_7%.xlsx');
data_O2_8 = readtable('DATA_8%.xlsx');
data_O2_9 = readtable('DATA_9%.xlsx');

tau_O2_6 = data_O2_6.Tau_value;
temperature_O2_6 = data_O2_6.Temperature;

tau_O2_7 = data_O2_7.Tau_value;
temperature_O2_7 = data_O2_7.Temperature;

tau_O2_8 = data_O2_8.Tau_value;
temperature_O2_8 = data_O2_8.Temperature;

tau_O2_9 = data_O2_9.Tau_value;
temperature_O2_9 = data_O2_9.Temperature;

calibrated_tau = tau_temp_calibration(tau, temperature);

% Apply the calibration function to each tau value
calibrated_tau_O2_6 = arrayfun(@(t, temp) tau_temp_calibration(t, temp), tau_O2_6, temperature_O2_6);
calibrated_tau_O2_7 = arrayfun(@(t, temp) tau_temp_calibration(t, temp), tau_O2_7, temperature_O2_7);
calibrated_tau_O2_8 = arrayfun(@(t, temp) tau_temp_calibration(t, temp), tau_O2_8, temperature_O2_8);
calibrated_tau_O2_9 = arrayfun(@(t, temp) tau_temp_calibration(t, temp), tau_O2_9, temperature_O2_9);

% 1. Calculate the Moving Average of the Calibrated Tau Value with a Window of 100
window_size = 100;
moving_avg_tau_O2_6 = movmean(calibrated_tau_O2_6, window_size);
moving_avg_tau_O2_7 = movmean(calibrated_tau_O2_7, window_size);
moving_avg_tau_O2_8 = movmean(calibrated_tau_O2_8, window_size);
moving_avg_tau_O2_9 = movmean(calibrated_tau_O2_9, window_size);

% 2. Plot the Moving Average of Calibrated Tau with Respect to Temperature
figure;
hold on;
% plot(temperature_O2_6, moving_avg_tau_O2_6, 'r-', 'LineWidth', 1.5, 'DisplayName', 'O2 = 6%');
plot(temperature_O2_7, moving_avg_tau_O2_7, 'r-', 'LineWidth', 2, 'DisplayName', 'O2 = 7%');
plot(temperature_O2_8, moving_avg_tau_O2_8, 'b-', 'LineWidth', 2, 'DisplayName', 'O2 = 8%');
% plot(temperature_O2_9, moving_avg_tau_O2_9, 'm-', 'LineWidth', 1.5, 'DisplayName', 'O2 = 9%');

xlabel('Temperature (°C)');
ylabel('Moving Average of Calibrated Tau (µs)');
title('Moving Average of Calibrated Tau with Respect to Temperature');
grid on;
legend show;

% 2. Plot the Moving Average of Calibrated Tau with Respect to Temperature
figure;
hold on;
% plot(temperature_O2_6, moving_avg_tau_O2_6, 'r.', 'LineWidth', 1.5, 'DisplayName', 'O2 = 6%');
plot(temperature_O2_7, calibrated_tau_O2_7 , 'r.', 'LineWidth', 2, 'DisplayName', 'O2 = 7%');
plot(temperature_O2_8, calibrated_tau_O2_8, 'b.', 'LineWidth', 2, 'DisplayName', 'O2 = 8%');
% plot(temperature_O2_9, moving_avg_tau_O2_9, 'm-', 'LineWidth', 1.5, 'DisplayName', 'O2 = 9%');

xlabel('Temperature (°C)');
ylabel('Moving Average of Calibrated Tau (µs)');
title('Moving Average of Calibrated Tau with Respect to Temperature');
grid on;
legend show;

% 3. Calculate the Calibrated Tau Value at 25 Degrees by Averaging from 24 to 26 Degrees
tau_25_O2_6 = mean(calibrated_tau_O2_6(temperature_O2_6 >= 24 & temperature_O2_6 <= 26));
tau_25_O2_7 = mean(calibrated_tau_O2_7(temperature_O2_7 >= 24 & temperature_O2_7 <= 26));
tau_25_O2_8 = mean(calibrated_tau_O2_8(temperature_O2_8 >= 24 & temperature_O2_8 <= 26));
tau_25_O2_9 = mean(calibrated_tau_O2_9(temperature_O2_9 >= 24 & temperature_O2_9 <= 26));

% Calculate the error with respect to the calibrated tau value at 24-26°C
tau_error_O2_6 = ((moving_avg_tau_O2_6 - tau_25_O2_6) / tau_25_O2_6) * 100;
tau_error_O2_7 = ((moving_avg_tau_O2_7 - tau_25_O2_7) / tau_25_O2_7) * 100;
tau_error_O2_8 = ((moving_avg_tau_O2_8 - tau_25_O2_8) / tau_25_O2_8) * 100;
tau_error_O2_9 = ((moving_avg_tau_O2_9 - tau_25_O2_9) / tau_25_O2_9) * 100;

% Plot the Error
figure;
hold on;
% plot(temperature_O2_6, tau_error_O2_6, 'r-', 'DisplayName', 'O2 = 6%');
plot(temperature_O2_7, tau_error_O2_7, 'r-', 'LineWidth', 2, 'DisplayName', 'O2 = 7%');
plot(temperature_O2_8, tau_error_O2_8, 'b-', 'LineWidth', 2, 'DisplayName', 'O2 = 8%');
% plot(temperature_O2_9, tau_error_O2_9, 'm-', 'DisplayName', 'O2 = 9%');

xlabel('Temperature (°C)');
ylabel('Error percentage in Calibrated Tau (%)');
title('Error of Calibrated Tau with Respect to Averaged Calibrated Tau at 24-26°C');
grid on;
legend show;

% 4. Apply the Least Squares Method to the Calibrated Tau Error
% Fit a linear model (y = a*x + b) to the calibrated tau error

X = [ones(length(temperature_O2_8), 1), temperature_O2_8]; % Design matrix for linear regression
coefficients = X \ tau_error_O2_8; % Least squares solution

% Predicted values using the linear fit
tau_error_fit = X * coefficients;

% Plot the fitted line against the original data
figure;
plot(temperature_O2_8, tau_error_O2_8, 'o'); % Original data
hold on;
plot(temperature_O2_8, tau_error_fit, '-r'); % Fitted line
xlabel('Temperature (°C)');
ylabel('Error in Calibrated Tau (µs)');
title('Least Squares Fit to Calibrated Tau Error');
legend('Tau Error', 'Least Squares Fit');
grid on;
hold off;

% Display the coefficients
a = coefficients(2); % Slope
b = coefficients(1); % Intercept
disp(['Fitted line: Tau Error = ' num2str(a) '*Temperature + ' num2str(b)]);



% Define the calibration function
function calibrated_tau = tau_temp_calibration(tau, temperature)
    % Initialize the scaling factor
    scaling_factor = 1.0;

    % Apply scaling factor based on temperature range
    if 20 < temperature <= 22.0
        scaling_factor = 0.95;
    elseif 22 < temperature <= 24.0
        scaling_factor = 1.0;
    elseif 24 < temperature <= 26.0
        scaling_factor = 1.02;
    elseif 26 < temperature <= 28.0
        scaling_factor = 1.04;
    elseif 28 < temperature <= 30.0
        scaling_factor = 1.06;
    elseif 30 < temperature <= 32.0
        scaling_factor = 1.08;
    elseif 32 < temperature <= 34.0
        scaling_factor = 1.10;
    elseif 34 < temperature <= 36.0
        scaling_factor = 1.12;
    elseif 36 < temperature <= 38.0
        scaling_factor = 1.14;
    elseif 38 < temperature <= 40.0
        scaling_factor = 1.16;
    elseif  40 < temperature <= 42.0
        scaling_factor = 1.18;
    elseif 42 < temperature <= 44.0
        scaling_factor = 1.2;
    end

    % Calculate the calibrated tau value
    calibrated_tau = tau * scaling_factor;
end