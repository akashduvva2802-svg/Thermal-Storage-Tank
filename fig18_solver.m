% fig18_solver.m
% Replicates Figure 18: Shifting the Hot Heat Exchanger vertically.

model = model_data();
tank = tank_data();

% Baseline parameters
model.Reh = 6.48e4; 
model.Rec = 1.73e4;

% Setup coil geometry for Figure 18
model.l_h = 1.0; % Hot coil is 1m long
model.l_c = 2.0; % Cold coil covers the entire tank (2m long)
model.x_c = 0.0; % Cold coil starts at the top (0m)

% Define the three configurations for the hot coil location
xh_configs = [0.1, 0.5, 0.9];
config_labels = {'x_h = 0.1m - 1.1m (Top)', 'x_h = 0.5m - 1.5m (Middle)', 'x_h = 0.9m - 1.9m (Bottom)'};

% Simulation settings
t_end = 10 * 3600; % 10 hours for steady state
dt_step = 60;      % 60-second steps for buoyancy mixing

Results_T_out = cell(length(xh_configs), 1);
node_vol = (pi*(tank.D^2)/4) * model.dx;

figure('Position', [100, 100, 800, 500]);
hold on; box on;

markers = {'d-', '*-', 'o-'};
h_sim = zeros(1, length(xh_configs));

for k = 1:length(xh_configs)
    disp(['Running simulation for hot coil at ', config_labels{k}]);
    
    model.x_h = xh_configs(k);
    
    % Initial conditions
    T_in = ones(model.N, 1) * model.T0;
    Th_in = ones(model.N, 1) * model.T0;
    Tc_in = ones(model.N, 1) * model.T0;
    T_current = [T_in; Th_in; Tc_in];
    
    t_current = 0;
    
    % Buoyancy Mixing Loop
    while t_current < t_end
        t_span = [t_current, t_current + dt_step];
        
        % Using the standard verified physics model (change_gamma)
        [~, T_step] = ode15s(@(t, T_state) change_gamma(t, T_state, model, tank), t_span, T_current);
        
        T_latest = T_step(end, :)';
        T_tank = T_latest(1:model.N); 
        
        % Buoyancy Node Mixing (Crucial for Fig 18's reduced stratification!)
        for sweep = 1:model.N 
            for i = 2:model.N
                if T_tank(i) > T_tank(i-1)
                    prop = water_properties(T_tank(i));
                    prop_1 = water_properties(T_tank(i-1));
                    m = node_vol * prop.rho;
                    m1 = node_vol * prop_1.rho;
                    
                    T_mixed = (m*T_tank(i) + m1*T_tank(i-1))/(m + m1);
                    T_tank(i) = T_mixed;
                    T_tank(i-1) = T_mixed;
                end
            end
        end
        
        T_latest(1:model.N) = T_tank;
        T_current = T_latest;
        t_current = t_current + dt_step;
    end
    
    % Extract steady state temperatures
    T_tank_steady = T_current(1:model.N)';
    T_star_steady = (T_tank_steady - 300) / 100;
    x_H = (0.5 : 1 : model.N) / model.N; 
    
    idx_markers = 1:round(model.N/10):model.N;
    
    h_sim(k) = plot(x_H, T_star_steady, markers{k}, 'LineWidth', 1.5, ...
         'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none', ...
         'Color', 'k');
end

xlabel('x/H', 'FontSize', 14);
ylabel('Dimensionless Temperature (T^*)', 'FontSize', 14);
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1.0, 'FontSize', 14);
axis([0 1 0 1.0]); 

lgd = legend(h_sim, config_labels);
set(lgd, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
title('Fig 18: Hot Coil Location & Buoyancy Effect', 'FontSize', 14);

hold off;
disp('Figure 18 generated successfully!');
