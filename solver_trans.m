t_span = 0:100:5400; % Time array up to 3 hours with 45s steps
model = model_data();
tank = tank_data();
T_in = ones(model.N, 1) * model.T0;
Th_in = ones(model.N, 1) * model.T0;
Tc_in = ones(model.N, 1) * model.T0;
T_initial = [T_in; Th_in; Tc_in];

% Solve ODE
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
[t_out, T_out] = ode15s(@(t, T_state) transient_ht(t, T_state, model, tank), t_span, T_initial, options);

% Extract raw temperatures
N = model.N;
T_raw  = T_out(:, 1:N);
Th_raw = T_out(:, N+1:2*N);
Tc_raw = T_out(:, 2*N+1:3*N);

% Dimensionless Temperatures
T_t0 = 300; 
delta_T_in = 100;

T_star  = (T_raw  - T_t0)  ./ delta_T_in;
Th_star = (Th_raw - T_t0)./delta_T_in;
Tc_out = T_out(:, 2*N + 1);
Tc_out_star = (Tc_out - T_t0) ./ delta_T_in;
Th_in_star = (function_Th_in(t_out, 1/1800) - T_t0) ./ delta_T_in;
T_avg_star = (sum(T_raw, 2) / N - T_t0) ./ delta_T_in;
X = t_out/3600;
Y1_ref = -42.9052.*X.^9 + 324.505.*X.^8 - 1022.08.*X.^7 + 1733.4.*X.^6 - 1707.13.*X.^5 + 978.974.*X.^4 - 308.137.*X.^3 + 44.9612.*X.^2 - 1.29812.*X - 0.0041152;
Y2_ref = -14.906.*X.^9 + 119.628.*X.^8 - 398.182.*X.^7 + 712.251.*X.^6 - 740.413.*X.^5 + 450.396.*X.^4 - 152.372.*X.^3 + 24.8225.*X.^2 - 1.04899.*X + 0.00124559;
Y3_ref = -9.8001.*X.^9 + 75.6742.*X.^8 - 243.578.*X.^7 + 422.871.*X.^6 - 427.568.*X.^5 + 253.043.*X.^4 - 83.0324.*X.^3 + 12.9656.*X.^2 - 0.459585.*X - 0.0026006;
Y4_ref = -13.2101.*X.^9 + 103.4.*X.^8 - 335.933.*X.^7 + 585.761.*X.^6 - 591.115.*X.^5 + 345.964.*X.^4 - 110.488.*X.^3 + 16.1488.*X.^2 - 0.344538.*X + 0.00613487;
Y5_ref = -16.899.*X.^9 + 136.135.*X.^8 - 455.345.*X.^7 + 819.454.*X.^6 - 858.344.*X.^5 + 527.423.*X.^4 - 181.135.*X.^3 + 30.3495.*X.^2 - 1.47297.*X - 0.00642475;

node_top = max(1, round(0.05 * N));
node_mid = max(1, round(0.55 * N));
node_bot = max(1, round(0.95 * N));

% Find indices for markers every 0.125 hours (which is 450 seconds)
% Since we used dt=45, 450s is every 10 indices
% However ode15s may change time steps, so we interpolate or just find closest
idx_markers = [];
target_times = 0:0.125:3;
for i = 1:length(target_times)
    [~, min_idx] = min(abs(X - target_times(i)));
    idx_markers = [idx_markers, min_idx];
end

% Plot exactly like Figure 7
figure('Position', [100, 100, 1200, 450]);
subplot(1, 2, 1);
hold on; box on;
% Plot lines and markers
p1 = plot(X, Th_in_star, 'k-', 'LineWidth', 1.2, 'Marker', '*', 'MarkerIndices', idx_markers, 'MarkerSize', 8);
p2 = plot(X, T_star(:, node_top), 'k-', 'LineWidth', 1.2, 'Marker', 'o', 'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none');
p3 = plot(X, T_star(:, node_mid), 'k-', 'LineWidth', 1.2, 'Marker', 's', 'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none');
p4 = plot(X, T_star(:, node_bot), 'k-', 'LineWidth', 1.2, 'Marker', 'd', 'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none');

% Plot the user polynomial as a dotted line without adding to legend
plot(X, Y1_ref, 'k:', 'LineWidth', 1.5);
plot(X, Y2_ref, 'k:', 'LineWidth', 1.5);
plot(X, Y3_ref, 'k:', 'LineWidth', 1.5);

xlabel('time (hours)', 'FontSize', 14);
ylabel('T^*', 'FontSize', 14);
set(gca, 'XTick', 0:1:3, 'YTick', 0:0.5:1.5, 'FontSize', 12);
axis([0 1.5 0 1.5]);

lgd = legend([p1, p2, p3, p4], {'Inlet hot water', 'T(x/H = 0.05)', 'T(x/H = 0.55)', 'T(x/H = 0.95)'});
set(lgd, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);

subplot(1, 2, 2);
hold on; box on;
p1_2 = plot(X, Th_in_star, 'k-', 'LineWidth', 1.2, 'Marker', 'o', 'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none');
p5 = plot(X, T_avg_star, 'k-', 'LineWidth', 1.2, 'Marker', 's', 'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none');
p6 = plot(X, Tc_out_star, 'k-', 'LineWidth', 1.2, 'Marker', 'd', 'MarkerIndices', idx_markers, 'MarkerSize', 8, 'MarkerFaceColor', 'none');

% Plot the user polynomial as a dotted line without adding to legend
plot(X, Y4_ref, 'k:', 'LineWidth', 1.5);
plot(X, Y5_ref, 'k:', 'LineWidth', 1.5);

xlabel('time (hours)', 'FontSize', 14);
ylabel('Dimensionless temperature', 'FontSize', 14);
set(gca, 'XTick', 0:1:3, 'YTick', 0:0.2:1.6, 'FontSize', 12);
axis([0 1.5 0 1.6]);

lgd2 = legend([p1_2, p5, p6], {'T^*_{h,in}', 'T^*_{avg}', 'T^*_{c,out}'});
set(lgd2, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);


R1_matrix = corrcoef(T_star(:, node_top), Y1_ref);
R1 = R1_matrix(1, 2);
R1_sq = R1^2;
MSE1 = mean((T_star(:, node_top) - Y1_ref).^2);
fprintf('x/H = 0.05 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R1, R1_sq, MSE1);
R2_matrix = corrcoef(T_star(:, node_mid), Y2_ref);
R2 = R2_matrix(1, 2);
R2_sq = R2^2;
MSE2 = mean((T_star(:, node_mid) - Y2_ref).^2);
fprintf('x/H = 0.55 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R2, R2_sq, MSE2);
R3_matrix = corrcoef(T_star(:, node_bot), Y3_ref);
R3 = R3_matrix(1, 2);
R3_sq = R3^2;
MSE3 = mean((T_star(:, node_bot) - Y3_ref).^2);
fprintf('x/H = 0.95 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R3, R3_sq, MSE3);

R4_matrix = corrcoef(T_avg_star, Y4_ref);
R4 = R4_matrix(1, 2);
R4_sq = R4^2;
MSE4 = mean((T_avg_star - Y4_ref).^2);
fprintf('T_avg_star : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R4, R4_sq, MSE4);

R5_matrix = corrcoef(Tc_out_star, Y5_ref);
R5 = R5_matrix(1, 2);
R5_sq = R5^2;
MSE5 = mean((Tc_out_star - Y5_ref).^2);
fprintf('Tc_out_star : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R5, R5_sq, MSE5);
