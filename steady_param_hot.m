% 1. Load the base data once
model = model_data();
tank = tank_data();

% 2. Define the array of Reynolds numbers you want to sweep (e.g., from Fig 10)
model.Reh_list = [2.16e3, 6.48e4, 1.6e5, 2.4e5];
Reh_array = model.Reh_list;

% Removed empirical tuning factors to preserve strictly rigorous theoretical modeling

% 3. Create Cell Arrays to store the results of each run
Results_T_out = cell(length(Reh_array), 1);
Results_t_out = cell(length(Reh_array), 1);

t_span = linspace(0, 3*3600, 200);
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
    T_out_run = Results_T_out{k};
    T_steady_all = T_out_run(end, :);
    T_tank_steady = T_steady_all(1:model.N);
    
    % Convert to Dimensionless Temperature
    T_star_steady = (T_tank_steady - 300) / 100;
    x_H = (0.5 : 1 : model.N) / model.N; 
    
    idx_markers = 1:round(model.N/10):model.N;
    
    % Plot the line for run 'k'
    h_sim(k) = plot(x_H, T_star_steady, ['k-', markers{k}], 'LineWidth', 1.2, ...
         'MarkerIndices', idx_markers, 'MarkerSize', 10, 'MarkerFaceColor', 'none');
         
    if k == 1
        Y_ref = 1122.56.*x_H.^9 - 4932.49.*x_H.^8 + 9167.89.*x_H.^7 - 9393.2.*x_H.^6 + 5795.8.*x_H.^5 - 2209.39.*x_H.^4 + 512.726.*x_H.^3 - 67.5101.*x_H.^2 + 3.36193.*x_H + 0.508186;
    elseif k == 2
        Y_ref = 1092.45.*x_H.^9 - 4830.81.*x_H.^8 + 9058.93.*x_H.^7 - 9396.16.*x_H.^6 + 5894.23.*x_H.^5 - 2295.17.*x_H.^4 + 546.754.*x_H.^3 - 74.901.*x_H.^2 + 4.46069.*x_H + 0.582742;
    elseif k == 3
        Y_ref = 990.913.*x_H.^9 - 4532.29.*x_H.^8 + 8808.86.*x_H.^7 - 9484.4.*x_H.^6 + 6180.31.*x_H.^5 - 2498.28.*x_H.^4 + 616.107.*x_H.^3 - 87.1848.*x_H.^2 + 5.8015.*x_H + 0.638736;
    elseif k == 4
        Y_ref = 154.471.*x_H.^9 - 710.429.*x_H.^8 + 1388.61.*x_H.^7 - 1505.*x_H.^6 + 989.877.*x_H.^5 - 406.433.*x_H.^4 + 103.099.*x_H.^3 - 15.4592.*x_H.^2 + 0.936066.*x_H + 0.776797;
    end
    
    plot(x_H, Y_ref, 'k:', 'LineWidth', 1.5);
    
    R_matrix = corrcoef(T_star_steady(:), Y_ref(:));
    R_val = R_matrix(1, 2);
    R_sq = R_val^2;
    MSE_val = mean((T_star_steady(:) - Y_ref(:)).^2);
    fprintf('Re_h = %e : Correlation (R) = %f, R^2 = %f, MSE = %e\n', Reh_array(k), R_val, R_sq, MSE_val);
end
% 3. Add your labels and legend AFTER the loop is finished
xlabel('x/H', 'FontSize', 14);
ylabel('T^*', 'FontSize', 14);
set(gca, 'XTick', 0:0.5:1, 'YTick', 0:0.2:1.6, 'FontSize', 14);
axis([0 1 0 1.7]);

% Add a legend to tell the 4 lines apart!
lgd = legend(h_sim, {['Re_h = 2.17 \times 10^3'], ...
              ['Re_h = 6.48 \times 10^4'], ...
              ['Re_h = 1.62 \times 10^5'], ...
              ['Re_h = 2.43 \times 10^5']});
set(lgd, 'Location', 'northeast', 'Box', 'off', 'FontSize', 14);
hold off;


