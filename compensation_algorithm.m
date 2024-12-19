%

%
% Load the data from Excel files
temp_23_data = readtable('ALL_O2_temp_room.xlsx');
temp_10_data = readtable('ALL_O2_temp10.xlsx');
temp_30_data = readtable('ALL_O2_temp30.xlsx');
temp_40_data = readtable('All_O2_temp40.xlsx');
% Load the data from Excel files

Tau_23 = temp_23_data.Tau_value_mean;
%Tau_10 = temp_10_data.Tau_value_mean;
Tau_30 = temp_30_data.Tau_value_mean;
Tau_40 = temp_40_data.Tau_value_mean;

O2mmHg_23 = temp_23_data.O2_mmHg; 
%O2mmHg_10 = temp_10_data.O2_mmHg;
O2mmHg_30 = temp_30_data.O2_mmHg;
O2mmHg_40 = temp_40_data.O2_mmHg;

% Applying stern-vlomer model

% Calculate Tau_value(O2_mmHg == 0) for each dataset
tau_value_0_23 = Tau_23(O2mmHg_23 == 0);
%tau_value_0_10 = Tau_10(O2mmHg_10 == 0);
tau_value_0_30 = Tau_30(O2mmHg_30 == 0);
tau_value_0_40 = Tau_40(O2mmHg_40 == 0);




% Calculate the ratio (Tau_value(O2_mmHg = 0) / Tau_value)-1 for each dataset
KSVO2_23 = (tau_value_0_23 ./ Tau_23)-1;
%KSVO2_10 = (tau_value_0_10 ./ Tau_10)-1;
KSVO2_30 = (tau_value_0_30 ./ Tau_30)-1;
KSVO2_40 = (tau_value_0_40 ./ Tau_40)-1;

% Plot the ratio with respect to O2_mmHg with fitted linear lines
figure;
hold on;
plot(O2mmHg_23, KSVO2_23 ,'o--b', 'Linewidth', 1.5, 'DisplayName', 'Temperature: 23°C');
%plot(O2mmHg_10, ratio_10 , '*-', 'Linewidth', 1.5, 'DisplayName', 'Temp 10°C');
plot(O2mmHg_30, KSVO2_30 , 'd--r', 'Linewidth', 1.5, 'DisplayName', 'Temperature: 30°C');
plot(O2mmHg_40, KSVO2_40 , '^--m', 'Linewidth', 1.5, 'DisplayName', 'Temperature: 40°C');
legend show;
xlabel('PtcO_2(mmHg)');
xticks(0:7.6:150) 
ylabel('(\tau_{0} / \tau)-1');
title('Stern-volmer correlation curves');
legend show;
grid on;



fit_type = fittype('KSV*x', 'independent', 'x', 'coefficients', 'KSV');

% 
% % Fit a linear model (poly1) to the ratio data
fit_linear_23 = fit(O2mmHg_23(1:10), KSVO2_23(1:10), fit_type);
%fit_linear_10 = fit(O2mmHg_10(1:10), KSVO2_10(1:10), fit_type);
fit_linear_30 = fit(O2mmHg_30(1:10), KSVO2_30(1:10), fit_type);
fit_linear_40 = fit(O2mmHg_40(1:10), KSVO2_40(1:10), fit_type);

% Calculate the residuals between the fitted curve and the actual data
%residuals_10 = ratio_10 - feval(fit_linear_10, O2mmHg_10);
residuals_23 = KSVO2_23 - feval(fit_linear_23, O2mmHg_23);
residuals_30 = KSVO2_30 - feval(fit_linear_30, O2mmHg_30);
residuals_40 = KSVO2_40 - feval(fit_linear_40, O2mmHg_40);

% Prevent Inf or NaN weights by adding a small constant
epsilon = 1e-6;
%weights_10 = 1 ./ (abs(residuals_10) + epsilon);
weights_23 = 1 ./ (abs(residuals_23) + epsilon);
weights_30 = 1 ./ (abs(residuals_30) + epsilon);
weights_40 = 1 ./ (abs(residuals_40) + epsilon);

% Apply the least squares method to minimize the sum of the squared residuals
% by refining the fit using the residuals
%lsq_fit_10 = fit(O2mmHg_10, ratio_10, fit_type, 'Weights', weights_10);
lsq_fit_23 = fit(O2mmHg_30, KSVO2_23, fit_type, 'Weights', weights_23);
lsq_fit_30 = fit(O2mmHg_30, KSVO2_30, fit_type, 'Weights', weights_30);
lsq_fit_40 = fit(O2mmHg_40, KSVO2_40, fit_type, 'Weights', weights_40);

fitted_values_23 = feval(lsq_fit_23, O2mmHg_23);
fitted_values_30 = feval(lsq_fit_30, O2mmHg_30);
fitted_values_40 = feval(lsq_fit_40, O2mmHg_40);

figure;
hold on;
plot(O2mmHg_23, fitted_values_23, 'o--b', 'Linewidth', 1.5, 'DisplayName', 'Temperature: 23°C ');
%plot(O2mmHg_10, normalized_ratio_10 , 'd-', 'Linewidth', 1.5,'DisplayName', 'Normalized Temp 10°C');
plot(O2mmHg_30, fitted_values_30, 'd--r', 'Linewidth', 1.5,'DisplayName', 'Temperature: 30°C');
plot(O2mmHg_40, fitted_values_40, '^--m','Linewidth', 1.5, 'DisplayName', 'Temperature: 40°C');
xlabel('PtcO_2 (mmHg)');
xticks(0:7.6:150) ;

ylabel('(\tau_{0} / \tau) -1 ');
title('stern-volmer correlation curves')
grid on;
legend show;


% Increase grid increments

% Calculate the scaling factors for the ratio curves
%scaling_factor_ratio_10 = mean(ratio_23) / mean(ratio_10);
scaling_factor_ratio_30 = mean(KSVO2_23) / mean(KSVO2_30);
scaling_factor_ratio_40 = mean(KSVO2_23) / mean(KSVO2_40);

% Normalize the ratio curves
%normalized_ratio_10 = ratio_10 * scaling_factor_ratio_10;
normalized_ratio_30 = KSVO2_30 * scaling_factor_ratio_30;
normalized_ratio_40 = KSVO2_40 * scaling_factor_ratio_40;


% Plot the normalized ratio curves with fitted linear lines
figure;
hold on;
plot(O2mmHg_23, KSVO2_23, 'o--b', 'Linewidth', 1.5, 'DisplayName', 'Temperature: 23°C ');
%plot(O2mmHg_10, normalized_ratio_10 , 'd-', 'Linewidth', 1.5,'DisplayName', 'Normalized Temp 10°C');
plot(O2mmHg_30, normalized_ratio_30, 'd--r', 'Linewidth', 1.5,'DisplayName', 'Normalized Temperature: 30°C');
plot(O2mmHg_40, normalized_ratio_40, '^--m','Linewidth', 1.5, 'DisplayName', 'Normalized Temperature: 40°C');
xlabel('PtcO_2 (mmHg)');
xticks(0:7.6:150) ;

ylabel('Normalized (\tau_{0} / \tau) -1 ');
title('Normalized the stern-volmer correlation curves')
grid on;
legend show;

% Create a string with the scaling factors
scaling_factors_ratio_text = sprintf('Avg Scaling Factor 30°C = %.2f\nAvg Scaling Factor 40°C = %.2f', ...
scaling_factor_ratio_30, scaling_factor_ratio_40);

% Display the scaling factors in a box 
annotation('textbox', [0.1, 0.1, 0.1, 0.1], 'String', scaling_factors_ratio_text, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

% **********************************************************************************
% Define maximum number of iterations and tolerance for convergence
tolerance = 0.01;
max_iterations = 10;

format long;
temp_30_data.Tau = unique (Tau_30);
temp_40_data.Tau = unique (Tau_40);
room_temp_data.Tau = unique (Tau_23);

% Initialize variables for iteration
%normalized_tau_10 = Tau_10;  % Initial guess (you can start with original Tau)
normalized_tau_30 = Tau_30;  % Same for temp 30
normalized_tau_40 = Tau_40;  % Same for temp 40

%prev_normalized_tau_10 = normalized_tau_10;
prev_normalized_tau_30 = normalized_tau_30;
prev_normalized_tau_40 = normalized_tau_40;

%scaling_factor_10 = 1;
scaling_factor_30 = 1;
scaling_factor_40 = 1;

for iter = 1:max_iterations
    % Step 1: Calculate scaling factors for each temperature
    %pre_scaling_factor_10 = Tau_23(1:11) ./ prev_normalized_tau_10;
    pre_scaling_factor_30 = Tau_23 ./ prev_normalized_tau_30;
    pre_scaling_factor_40 = Tau_23 ./ prev_normalized_tau_40;

    %scaling_factor_10 = pre_scaling_factor_10 .* scaling_factor_10;
    scaling_factor_30 = pre_scaling_factor_30 .* scaling_factor_30;
    scaling_factor_40 = pre_scaling_factor_40 .* scaling_factor_40;

    % Step 2: Update the normalized tau values
    %normalized_tau_10 = Tau_10 .* mean(scaling_factor_10);
    normalized_tau_30 = Tau_30 .* mean(scaling_factor_30);
    normalized_tau_40 = Tau_40 .* mean(scaling_factor_40);

    % Step 3: Check for convergence (difference between current and previous)
    %diff_10 = min(abs(normalized_tau_10 - Tau_23(1:11)));
    diff_30 = min(abs(normalized_tau_30 - Tau_23));
    diff_40 = min(abs(normalized_tau_40 - Tau_23));

    % Update previous values for the next iteration
    %prev_normalized_tau_10 = normalized_tau_10;
    prev_normalized_tau_30 = normalized_tau_30;
    prev_normalized_tau_40 = normalized_tau_40;



    % Check if the changes are below the tolerance
    if  diff_30 < tolerance && diff_40 < tolerance
        disp(['Convergence reached after ' num2str(iter) ' iterations.']);
        break;
    end
end

% If maximum iterations reached
if iter == max_iterations
    disp('Maximum number of iterations reached without full convergence.');
end

% Add normalized tau values to the respective data tables
%temp_10_data.Normalized_Tau = normalized_tau_10;
temp_30_data.Normalized_Tau = normalized_tau_30;
temp_40_data.Normalized_Tau = normalized_tau_40;

% Plot the normalized tau values
figure;
hold on;
plot(O2mmHg_23, Tau_23, 'o--', 'Linewidth', 1.5, 'DisplayName', 'Temperature: 23°C ');
%plot(O2mmHg_10, normalized_tau_10, 'd', 'Linewidth', 1.5, 'DisplayName', 'Normalized Temp 10°C');
plot(O2mmHg_30, normalized_tau_30, 'd--', 'Linewidth', 1.5, 'DisplayName', 'Normalized Temperature: 30°C');
plot(O2mmHg_40, normalized_tau_40, '^--', 'Linewidth', 1.5, 'DisplayName', 'Normalized Temperature: 40°C');
legend show;
xlabel('O2 partial pressure (mmHg)');
xticks(0:7.6:150) 

ylabel('\tau Value (\mu s)');
title('Iterative Normalization of \tau Values');
grid on;


% scaling_factor_10 = Tau_23(1:11) ./ Tau_10;   
% scaling_factor_30 = Tau_23 ./ Tau_30;
% scaling_factor_40 = Tau_23 ./ Tau_40;
% 
% % Normalize the tau values
% normalized_tau_10 = Tau_10 .* mean(scaling_factor_10);
% normalized_tau_30 = Tau_30 .* mean(scaling_factor_30);
% normalized_tau_40 = Tau_40 .* mean(scaling_factor_40);
% 
% % Add normalized tau values to the respective data tables
% temp_10_data.Normalized_Tau = normalized_tau_10;
% temp_30_data.Normalized_Tau = normalized_tau_30;
% temp_40_data.Normalized_Tau = normalized_tau_40;

% Fit exponential model (exp2) to the data
fit_exp2_23 = fit(O2mmHg_23, Tau_23, 'exp2');
fit_exp2_30 = fit(O2mmHg_30, Tau_30, 'exp2');
fit_exp2_40 = fit(O2mmHg_40, Tau_40, 'exp2');

% Fit exponential model (exp2) to the data after normalization
fit_exp2_30_norm = fit(O2mmHg_30, normalized_tau_30, 'exp2');
fit_exp2_40_norm = fit(O2mmHg_40, normalized_tau_40, 'exp2');

% Plot the tau values before normalization with fitted models
figure;
hold on;
legend show;

plot(O2mmHg_23, Tau_23, 'o', 'Linewidth', 1.5, 'DisplayName', 'Temperature: 23°C');
%plot(O2mmHg_10, Tau_10, 'd','Linewidth', 1.5, 'DisplayName', 'Temp 10°C');
plot(O2mmHg_30, Tau_30, 'd','Linewidth', 1.5, 'DisplayName', 'Temperature: 30°C');
plot(O2mmHg_40, Tau_40, '^m','Linewidth', 1.5, 'DisplayName', 'Temperature: 40°C');

plot(fit_exp2_23, 'b--');
plot(fit_exp2_30, 'r--');
plot(fit_exp2_40, 'm--');
legend off;

xlabel('O2 partial pressure (mmHg)');
ylabel('\tau (\mu s)');
title('\tau Values Before Normalization with Fitted Models');
xticks(0:7.6:150) 
grid on;

% Plot the tau values after normalization with fitted models
figure;
hold on;
plot(O2mmHg_23, Tau_23, 'o--b', 'Linewidth', 1.5,'DisplayName', 'Temperature: 23°C');
%plot(O2mmHg_10, temp_10_data.Normalized_Tau,  'd-', 'Linewidth', 1.5, 'DisplayName', 'Normalized Temp 10°C');
plot(O2mmHg_30, temp_30_data.Normalized_Tau,  'd--r', 'Linewidth', 1.5, 'DisplayName', 'Normalized Temp 30°C');
plot(O2mmHg_40, temp_40_data.Normalized_Tau, '^--m','Linewidth', 1.5, 'DisplayName', 'Normalized Temp 40°C');
legend show;
% plot(fit_exp2_23,'b--');
% plot(fit_exp2_30_norm, 'y--');
% plot(fit_exp2_40_norm, 'm--');
xlabel('PtcO_2 (mmHg)');
xticks(0:7.6:150) ;
ylabel('\tau  (\mu s)');
title('\tau After Normalization');
grid on;


% Create a string with the scaling factors
scaling_factors_text = sprintf('Avg Scaling Factor 30°C = %.2f\nAvg Scaling Factor 40°C = %.2f', ...
    mean(scaling_factor_30), mean(scaling_factor_40));

% Display the scaling factors in a box 
annotation('textbox', [0.5, 0.5, 0.1, 0.1], 'String', scaling_factors_text, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

% text(1, max(Tau_23) - 1, ['Avg Scaling Factor 30°C = ', num2str(mean(scaling_factor_30))], 'FontSize', 12, 'Color', 'k');
% text(1, max(Tau_23) - 2, ['Avg Scaling Factor 40°C = ', num2str(mean(scaling_factor_40))], 'FontSize', 12, 'Color', 'k');

%**************************************************************************

KSV_room = (KSVO2_23(O2mmHg_23 ~=0) - 1) ./ O2mmHg_23(2:end) ;
%KSV_10 = ((ratio_10(O2mmHg_10 ~=0) - 1) ./ O2mmHg_10(2:end)) ;
KSV_30 = ((KSVO2_30(O2mmHg_30 ~=0) - 1) ./ O2mmHg_30(2:end)) ;
KSV_40 = ((KSVO2_40(O2mmHg_40 ~=0) - 1) ./ O2mmHg_40(2:end)) ;


%KSV_10_avg = mean (KSV_10);
KSV_room_avg = mean (KSV_room);
KSV_30_avg = mean (KSV_30);
KSV_40_avg = mean (KSV_40);


%KSV_scaling_factor_10 = mean (KSV_room (1:10) ./ KSV_10);
KSV_scaling_factor_30 = mean (KSV_room ./ KSV_30);
KSV_scaling_factor_40 = mean (KSV_room ./ KSV_40);

%normalized_KSV_10 = KSV_scaling_factor_10 .* KSV_10;
normalized_KSV_30 = KSV_scaling_factor_30 .* KSV_30;
normalized_KSV_40 = KSV_scaling_factor_40 .* KSV_40;

%calibrated_ratio_10 (O2mmHg_10 ~= 0) = 1 + (normalized_KSV_10 .* O2mmHg_10 (2:end));
calibrated_ratio_30 (O2mmHg_30 ~= 0) = 1 + (normalized_KSV_30 .* O2mmHg_30(2:end));
calibrated_ratio_40 (O2mmHg_40 ~= 0) = 1 + (normalized_KSV_40 .* O2mmHg_40(2:end));
%calibrated_ratio_10 (O2mmHg_10 == 0) = 1;
calibrated_ratio_30 (O2mmHg_10 == 0) = 1;
calibrated_ratio_40 (O2mmHg_10 == 0) = 1;
figure;
hold on;
plot(O2mmHg_23, KSVO2_23 - 1, 'o-', 'Linewidth', 1.5,'DisplayName', 'Temperature: 23°C');
%plot(O2mmHg_10,calibrated_ratio_10 - 1, '*-', 'Linewidth', 1.5, 'DisplayName', 'Normalized Temp 10°C')
plot(O2mmHg_30, calibrated_ratio_30 - 1,  'd-', 'Linewidth', 1.5, 'DisplayName', 'Normalized Temp 30°C');
plot(O2mmHg_40, calibrated_ratio_40 - 1, '^-','Linewidth', 1.5, 'DisplayName', 'Normalized Temp 40°C');
legend show;
xlabel('O2 (mmHg)');
ylabel('(Tau_{0} / \tau)-1');
title(' Normalizing the Tau_{0} / \tau based on the Stern-Volmer Factor');
grid on;


%************************************************************************
% Extract unique O2 percentages
unique_O2 = unique(O2mmHg_23);

% Extract temperatures
temperatures_C = [23, 30, 40];

% Initialize a figure
figure;
hold on;

% Define a color map
colors = lines(length(unique_O2));

% Loop through each unique O2 percentage
for i = 1:length(unique_O2)
    O2_mmHg = unique_O2(i);
    
    % Filter data for the current O2 percentage
    Tau_10 = data10.Tau_value_mean(data10.O2_mmHg == O2_mmHg);
    Tau_23 = data22.Tau_value_mean(data22.O2_mmHg == O2_mmHg); 
    Tau_30 = data30.Tau_value_mean(data30.O2_mmHg == O2_mmHg); 
    Tau_40 = data40.Tau_value_mean(data40.O2_mmHg == O2_mmHg); 
    
    % Combine tau values for plotting
    tau_values = [Tau_23; Tau_30; Tau_40];
   
     % Fit a polynomial to the tau values vs temperature
    fit_result = fit(temperatures_C', tau_values, 'poly1'); % Linear fit 
    
    % Display the fit equation
    coeffs = coeffvalues(fit_result);
    fprintf('Fit equation for O2 = %d mmHg: tau = %.4f * T + %.4f\n', O2_mmHg, coeffs(1), coeffs(2));
    % Plot tau values against temperatures
    plot(temperatures_C, tau_values, '-o', 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', sprintf('O2 = %d mmHg', unique_O2(i)));
grid on;
legend('show');
xlabel('Temperature (°C)');
ylabel('Tau Value (us)');
title('Tau Value vs Temperature for Different O2 Percentages');
end
% Add labels and legend


hold off;





