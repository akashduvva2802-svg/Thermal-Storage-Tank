% Script to generate the plot for Figure 16 (Cold Heat Exchanger Location)
% This script sweeps different vertical locations of the 1m cold coil

% 1. Load the base data
model = model_data();
tank = tank_data();

% Set baseline Reynolds numbers for the parametric study (adjust if paper specifies differently)
model.Reh = 6.48e4; 
model.Rec = 1.73e4;

% Setup coil geometry for Figure 16
model.l_h = 2.0; % Hot coil is 2m long (covers entire tank)
model.x_h = 0.0; % Hot coil starts at the top (0m)
model.l_c = 1.0; % Cold coil is 1m long

% Define the three configurations for the cold coil location
% The paper evaluates top, middle, and bottom locations
xc_configs = [0.1, 0.5, 0.9];
config_labels = {'y = 0.1m - 1.1m (Top)', 'y = 0.5m - 1.5m (Middle)', 'y = 0.9m - 1.9m (Bottom)'};

% Create Cell Arrays to store the results
Results_T_out = cell(length(xc_configs), 1);
Results_t_out = cell(length(xc_configs), 1);

t_span = linspace(0, 5*3600, 200); % 10 hour simulation to ensure steady state

for k = 1:length(xc_configs)
    disp(['Running simulation for cold coil at ', config_labels{k}]);
    
    % Update the cold coil start location
    model.x_c = xc_configs(k);
    
    % Setup initial conditions (uniform at T0)
    T_in = ones(model.N, 1) * model.T0;
    Th_in = ones(model.N, 1) * model.T0;
    Tc_in = ones(model.N, 1) * model.T0;
    T_initial = [T_in; Th_in; Tc_in];
    
    % Run the ODE solver using the change_gamma model
    [t_out, T_out] = ode15s(@(t, T_state) change_gamma(t, T_state, model, tank), t_span, T_initial);
    
    % Save the data
    Results_t_out{k} = t_out;
    Results_T_out{k} = T_out;
end

% ---------------------------------------------------------
% Plotting the Results
% ---------------------------------------------------------
figure('Position', [100, 100, 600, 500]);
hold on; box on;

% The user's Curve 1 corresponds to Top [0.1, 1.1]
% The user's Curve 2 corresponds to Middle [0.5, 1.5]
% The user's Curve 3 corresponds to Bottom [0.9, 1.9]

% Curve 1 (User specifies this is now Bottom)
ref_T_bot = [0.9261, 0.9079, 0.8898, 0.8669, 0.7222, 0.6240, 0.5714, 0.5247, 0.4871, 0.4851];
% Curve 2 (User specifies this is Middle)
ref_T_mid = [0.9253, 0.9081, 0.7549, 0.6556, 0.6053, 0.5657, 0.5651, 0.5651, 0.5654, 0.5665];
% Curve 3 (User specifies this is now Top)
ref_T_top = [0.7763, 0.6771, 0.6288, 0.6033, 0.6047, 0.6058, 0.6051, 0.6063, 0.6074, 0.6068];

ref_T_all = {ref_T_top, ref_T_mid, ref_T_bot};

markers = {'d-', '*-', 'o-'};

h_sim = zeros(1, length(xc_configs));

for k = 1:length(xc_configs)
    % Extract the steady state data
    T_out_run = Results_T_out{k};
    T_steady_all = T_out_run(end, :);
    T_tank_steady = T_steady_all(1:model.N);
    
    % Convert to Dimensionless Temperature
    T_star_steady = (T_tank_steady - 300) / 100;
    x_H = (0.5 : 1 : model.N) / model.N; 
    
    idx_markers = 1:round(model.N/10):model.N;
    
    % Plot the simulation line
    h_sim(k) = plot(x_H, T_star_steady, markers{k}, 'LineWidth', 1.5, ...
         'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none', ...
         'Color', 'k');
         
    % Plot the reference data from the paper
    plot(x_H, ref_T_all{k}, [markers{k}(1), ':'], 'LineWidth', 1.5, 'Color', 'k', 'MarkerSize', 8);
    
    % Calculate exact R, R^2, and MSE
    Y_ref = ref_T_all{k};
    R_matrix = corrcoef(T_star_steady(:), Y_ref(:));
    R_val = R_matrix(1, 2);
    R_sq = R_val^2;
    MSE_val = mean((T_star_steady(:) - Y_ref(:)).^2);
    
    fprintf('%s : Correlation (R) = %f, R^2 = %f, MSE = %e\n', config_labels{k}, R_val, R_sq, MSE_val);
end

xlabel('x/H', 'FontSize', 14);
ylabel('T^*', 'FontSize', 14);
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1.0, 'FontSize', 14);
axis([0 1 0 1.0]); % Adjust Y-axis if necessary

% Add legend
lgd = legend(h_sim, config_labels);
set(lgd, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
title('Fig 16: Cold Coil Location Parametric Study', 'FontSize', 14);

hold off;
disp('Figure 16 with statistical correlations generated successfully!');
