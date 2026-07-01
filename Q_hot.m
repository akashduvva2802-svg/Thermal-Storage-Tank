% Q_hot.m
% Plot heat transferred from hot water vs time under different hot water Re

% 1. Load the base data once
model = model_data();
tank = tank_data();

% 2. Define the array of Reynolds numbers from the request
Reh_array = [2.17e3, 6.48e4, 1.62e5, 2.43e5];

% 3. Create Cell Arrays to store the results of each run
Results_T_out = cell(length(Reh_array), 1);
Results_t_out = cell(length(Reh_array), 1);

% Define time span: 3 hours
t_span = linspace(0, 3*3600, 10000); 

% Use Th_in = 400 for this specific plot
model.Th_in = 400;

for k = 1:length(Reh_array)
    disp(['Running simulation for Reh = ', num2str(Reh_array(k))]);
    
    % Override the hot water Reynolds number
    model.Reh = Reh_array(k);
    
    % Setup initial conditions (cold start, everything at T0)
    T_in = ones(model.N, 1) * model.T0;
    Th_in = ones(model.N, 1) * model.T0;
    Tc_in = ones(model.N, 1) * model.T0;
    T_initial = [T_in; Th_in; Tc_in];
    
    % Run the ODE solver
    [t_out, T_out] = ode15s(@(t, T_state) heat_trans_model(t, T_state, model, tank), t_span, T_initial);
    
    % Save results
    Results_t_out{k} = t_out;
    Results_T_out{k} = T_out;
end

% Plotting
figure('Position', [100, 100, 600, 500]);
hold on;
markers = {'d', '*', 'o', 's'};

for k = 1:length(Reh_array)
    t_out = Results_t_out{k};
    T_out_run = Results_T_out{k};
    
    % Extract hot water outlet temperature (node N of hot water, which is index 2*N)
    Th_out = T_out_run(:, 2*model.N);
    
    % Recalculate mass flow rate for this specific Reh
    model.Reh = Reh_array(k);
    prop_in_h = water_properties(model.Th_in);
    mh = model.Reh * (pi * tank.di * prop_in_h.mu) / 4;
    
    % Calculate Q in kW by integrating heat transfer from coil to tank
    Q_hot_kW = calc_Q_hot(T_out_run, model, tank);
    
    t_hours = t_out / 3600;
    
    % Plot with specific markers and color to match the requested graph
    plot(t_hours, Q_hot_kW, '-', 'Marker', markers{k}, 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerIndices', round(linspace(1, length(t_hours), 15)), 'Color', 'k');
end

% Formatting the plot to match the image
xlabel('time (hours)', 'FontSize', 12);
ylabel('Heat transferred from hot water (KW)', 'FontSize', 12);

% Legend entries with LaTeX formatting for scientific notation
legend_entries = {['Re_h = 2.17 \times 10^4'], ...
                  ['Re_h = 6.48 \times 10^4'], ...
                  ['Re_h = 1.62 \times 10^5'], ...
                  ['Re_h = 2.43 \times 10^5']};
legend(legend_entries, 'Location', 'northeast', 'FontSize', 12);

% Ensure limits and appearance are correct
xlim([0 3]);
ylim([0 300]); % As shown in the image
box on;
hold off;