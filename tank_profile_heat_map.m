% tank_profile_heat_map.m
% Generates colorful stratification heat maps for the 3 cold coil locations

model = model_data();
tank = tank_data();

% Set baseline Reynolds numbers as per user's original un-tuned base
model.Reh = 6.48e4; 
model.Rec = 1.73e4;

% Setup coil geometry for the study
model.l_h = 2.0; % Hot coil is 2m long (covers entire tank)
model.x_h = 0.0; % Hot coil starts at the top (0m)
model.l_c = 1.0; % Cold coil is 1m long

% Define the three configurations for the cold coil location
xc_configs = [0.1, 0.5, 0.9];
config_labels = {'Top Coil (0.1m - 1.1m)', 'Middle Coil (0.5m - 1.5m)', 'Bottom Coil (0.9m - 1.9m)'};

t_span = linspace(0, 10*3600, 200); % 10 hour simulation to reach steady state

% Store steady state data for the second plot
T_star_steady_all = cell(length(xc_configs), 1);

% Prepare wide figure for transient 3 subplots
fig1 = figure('Position', [100, 100, 1500, 500]);

for k = 1:length(xc_configs)
    disp(['Running simulation for ', config_labels{k}]);
    
    % Update the cold coil start location
    model.x_c = xc_configs(k);
    
    % Setup initial conditions (uniform at T0)
    T_in = ones(model.N, 1) * model.T0;
    Th_in = ones(model.N, 1) * model.T0;
    Tc_in = ones(model.N, 1) * model.T0;
    T_initial = [T_in; Th_in; Tc_in];
    
    % Run the ODE solver
    [t_out, T_out] = ode15s(@(t, T_state) change_gamma(t, T_state, model, tank), t_span, T_initial);
    
    % Extract Tank temperatures (Node 1 to N)
    T_tank = T_out(:, 1:model.N);
    
    % Convert to Dimensionless Temperature T*
    T_star = (T_tank - 300) / 100;
    
    % Save steady state for later
    T_star_steady_all{k} = T_star(end, :);
    
    % Setup meshgrid for contourf
    % Time in hours on X axis
    time_hrs = t_out / 3600;
    % Depth (dimensionless x/H) on Y axis. Node centers.
    depth_xH = (0.5 : 1 : model.N) / model.N; 
    
    [X, Y] = meshgrid(time_hrs, depth_xH);
    Z = T_star'; % Transpose because Z rows correspond to Y, columns to X
    
    % Create subplot
    figure(fig1);
    subplot(1, 3, k);
    contourf(X, Y, Z, 60, 'LineStyle', 'none');
    
    % Apply colorful colormap (jet is standard for thermal heat maps)
    colormap('jet');
    
    % Add colorbar and enforce fixed scale [0, 1] for direct visual comparison
    cb = colorbar;
    ylabel(cb, 'Dimensionless Temperature (T^*)', 'FontSize', 12);
    caxis([0 1]); 
    
    % Invert Y-axis so the top of the tank is at the top of the plot
    set(gca, 'YDir', 'reverse');
    
    % Draw dashed white lines to visually show exactly where the cold coil is located
    hold on;
    coil_start_y = xc_configs(k) / tank.H;
    coil_end_y = (xc_configs(k) + model.l_c) / tank.H;
    yline(coil_start_y, 'w--', 'LineWidth', 2);
    yline(coil_end_y, 'w--', 'LineWidth', 2);
    hold off;
    
    xlabel('Time (Hours)', 'FontSize', 14);
    ylabel('Depth (x/H)', 'FontSize', 14);
    title(config_labels{k}, 'FontSize', 16);
    
    % Minor formatting
    set(gca, 'FontSize', 12);
end

% Add a master title to figure 1
sgtitle('Transient Thermal Stratification Evolution', 'FontSize', 18, 'FontWeight', 'bold');


% --- Add Steady State Tank Slices ---
fig2 = figure('Position', [150, 150, 1000, 600]);

for k = 1:length(xc_configs)
    T_steady = T_star_steady_all{k}';
    
    % Create a 2D grid to represent the physical tank (Width vs Depth)
    % Width is arbitrary [0, 1] just for visualization
    width_x = linspace(0, 1, 20); 
    depth_y = (0.5 : 1 : model.N) / model.N;
    
    [X_ss, Y_ss] = meshgrid(width_x, depth_y);
    Z_ss = repmat(T_steady, 1, length(width_x));
    
    figure(fig2);
    subplot(1, 3, k);
    contourf(X_ss, Y_ss, Z_ss, 60, 'LineStyle', 'none');
    colormap('jet');
    caxis([0 1]);
    
    set(gca, 'YDir', 'reverse');
    
    % Draw coil locations
    hold on;
    coil_start_y = xc_configs(k) / tank.H;
    coil_end_y = (xc_configs(k) + model.l_c) / tank.H;
    
    % Draw lines to show the coil limits
    yline(coil_start_y, 'w--', 'LineWidth', 2);
    yline(coil_end_y, 'w--', 'LineWidth', 2);
    hold off;
    
    % Hide X-axis ticks because it's just a conceptual width
    set(gca, 'XTick', []);
    ylabel('Depth (x/H)', 'FontSize', 14);
    title(config_labels{k}, 'FontSize', 16);
    
    % Add colorbar to the last plot only
    if k == 3
        cb = colorbar;
        ylabel(cb, 'Dimensionless Temperature (T^*)', 'FontSize', 12);
    end
end

sgtitle('Steady-State Physical Tank Stratification Profiles', 'FontSize', 18, 'FontWeight', 'bold');

% Save the generated figures
saveas(fig1, 'tank_profile_heat_map_transient.png');
saveas(fig2, 'tank_profile_heat_map_steady.png');

disp('Both transient and steady-state heat maps generated and saved successfully!');

