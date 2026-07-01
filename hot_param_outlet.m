% 1. Load the base data once
model = model_data();
tank = tank_data();
N = model.N;
t_span = linspace(0, 3*3600, 200);
% 2. Define the array of Reynolds numbers you want to sweep (e.g., from Fig 10)
Reh_array = [2.17e3, 6.48e4, 1.62e5, 2.43e5];

% 3. Create Cell Arrays to store the results of each run
Results_T_out = cell(length(Reh_array), 1);
Results_t_out = cell(length(Reh_array), 1);

for k = 1:length(Reh_array)
    disp(['Running simulation for Reh = ', num2str(Reh_array(k))]);
    
    % --- THE MAGIC STEP: Override the Reynolds number ---
    model.Reh = Reh_array(k);
    
    % Setup initial conditions
    T_in = ones(model.N, 1) * model.T0;
    Th_in = ones(model.N, 1) * model.T0;
    Tc_in = ones(model.N, 1) * model.T0;
    T_initial = [T_in; Th_in; Tc_in];
    
    % Run the ODE solver. 
    % We use an anonymous function @(t, T_state) to pass our updated 'model' struct in!
    [t_out, T_out] = ode15s(@(t, T_state) heat_trans_model(t, T_state, model, tank), t_span, T_initial);
    
    % Save the data for this specific Reynolds number into the cell array
    Results_t_out{k} = t_out;
    Results_T_out{k} = T_out;
end

% 1. Create the figure window just ONCE
figure('Position', [100, 100, 600, 500]);
hold on; box on;

markers = {'d', '*', 'o', 's'};

% 2. Loop through all 4 of your saved results to plot them
h_sim = zeros(1, length(Reh_array));
for k = 1:length(Reh_array)
    
    % Extract the steady state data for run 'k'
    t_run = Results_t_out{k};
    T_out_run = Results_T_out{k};
    T_steady_all = T_out_run(:,2*N + 1);
    
    % Convert to Dimensionless Temperature
    T_star_steady = (T_steady_all - 300) / 100;
    
    % Find indices for markers every ~0.25 hours
    X = t_run / 3600;
    idx_markers = [];
    target_times = 0:0.25:3;
    for i = 1:length(target_times)
        [~, min_idx] = min(abs(X - target_times(i)));
        idx_markers = [idx_markers, min_idx];
    end
    
    if k == 1
        % Eq 4: lowest curve for lowest Re_h
        Y_ref = -0.00405325.*X.^5 + 0.0332407.*X.^4 - 0.095477.*X.^3 + 0.0855467.*X.^2 + 0.127931.*X - 0.0014563;
    elseif k == 2
        % Eq 3: 3rd highest curve for 2nd lowest Re_h
        Y_ref = -0.00373826.*X.^5 + 0.0315918.*X.^4 - 0.0889024.*X.^3 + 0.0462634.*X.^2 + 0.254723.*X - 6.73793e-05;
    elseif k == 3
        % Eq 2: 2nd highest curve for 2nd highest Re_h
        Y_ref = -0.00841815.*X.^5 + 0.0664149.*X.^4 - 0.165463.*X.^3 + 0.0448385.*X.^2 + 0.43763.*X - 0.000140872;
    elseif k == 4
        % Eq 1: highest curve for highest Re_h
        Y_ref = 0.000175774.*X.^4 + 0.025082.*X.^3 - 0.208939.*X.^2 + 0.605908.*X - 0.00211363;
    end
    
    % Plot the line for run 'k'
    h_sim(k) = plot(X, T_star_steady, ['k-', markers{k}], 'LineWidth', 1.2, ...
         'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none');
         
    plot(X, Y_ref, 'k:', 'LineWidth', 1.5);
    
    R_matrix = corrcoef(T_star_steady, Y_ref);
    R_val = R_matrix(1, 2);
    R_sq = R_val^2;
    MSE_val = mean((T_star_steady - Y_ref).^2);
    fprintf('Re_h = %e : Correlation (R) = %f, R^2 = %f, MSE = %e\n', Reh_array(k), R_val, R_sq, MSE_val);
end
% 3. Add your labels and legend AFTER the loop is finished
xlabel('time (hours)', 'FontSize', 14);
ylabel('Dimensionless cold water outlet temp (T^*)', 'FontSize', 14);
set(gca, 'XTick', 0:1:3, 'YTick', 0:0.2:1.6, 'FontSize', 12);
axis([0 3 0 1.7]);

% Add a legend to tell the 4 lines apart!
lgd = legend(h_sim, {['Re_h = 2.17 \times 10^3'], ...
              ['Re_h = 6.48 \times 10^4'], ...
              ['Re_h = 1.62 \times 10^5'], ...
              ['Re_h = 2.43 \times 10^5']});
set(lgd, 'Location', 'northeast', 'Box', 'off', 'FontSize', 14);
hold off;


