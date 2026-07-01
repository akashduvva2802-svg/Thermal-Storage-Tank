t_span = linspace(0, 3*3600, 10000);
model = model_data();
tank = tank_data();
T_in = ones(model.N, 1) * model.T0;
Th_in = ones(model.N, 1) * model.T0;
Tc_in = ones(model.N, 1) * model.T0;
T_initial = [T_in; Th_in; Tc_in];

[t_out, T_out] = ode15s(@(t, T_state) heat_trans_model(t, T_state, model, tank), t_span, T_initial);

% ==============================================================
% 1. Extract raw temperatures (Time x Nodes)
% ==============================================================
N = model.N;
T_raw  = T_out(:, 1:N);
Th_raw = T_out(:, N+1:2*N);
Tc_raw = T_out(:, 2*N+1:3*N);

% ==============================================================
% 2. Calculate Dimensionless Temperatures (Eq 27)
% ==============================================================
% Get the t=0 temperatures (the very first row of the matrix)
T_t0  = T_raw(1, :);
Th_t0 = Th_raw(1, :);
Tc_t0 = Tc_raw(1, :);

% Temperature denominator from Eq 27
delta_T_in = model.Th_in - model.Tc_in;

% Calculate T*, Th*, and Tc* for the entire matrix instantly
T_star  = (T_raw  - T_t0)  ./ delta_T_in;
Th_star = (Th_raw - Th_t0) ./ delta_T_in;
Tc_star = (Tc_raw - Tc_t0) ./ delta_T_in;

% ==============================================================
% 3. Plot the Dimensionless Results (Figs 4, 5, 6)
% ==============================================================
% Find the closest nodes to the paper's x/H = 0.05, 0.55, 0.95
node_top = max(1, round(0.05 * N));
node_mid = max(1, round(0.55 * N));
node_bot = max(1, round(0.95 * N));

t_hours = t_out / 3600;

figure('Position', [100, 100, 1200, 400]); % Make a wide window for 3 plots

% Plot T* (Stored Water - Fig 4)
subplot(1, 3, 1);
% Simulation plots (solid black lines)
plot(t_hours, T_star(:, node_top), 'k-', 'LineWidth', 1.5); hold on;
plot(t_hours, T_star(:, node_mid), 'k--', 'LineWidth', 1.5);
plot(t_hours, T_star(:, node_bot), 'k-.', 'LineWidth', 1.5);

% Reference polynomial plots
X = t_hours;
Y1_ref = -6.89645e-05*X.^9 + 0.00158494*X.^8 - 0.0152568*X.^7 + 0.0791985*X.^6 - 0.234743*X.^5 + 0.367097*X.^4 - 0.133142*X.^3 - 0.549789*X.^2 + 1.07224*X + 0.00770142;
Y2_ref = -2.58011e-06*X.^9 + 8.27623e-05*X.^8 - 0.00106777*X.^7 + 0.00739714*X.^6 - 0.0304093*X.^5 + 0.0744141*X.^4 - 0.085256*X.^3 - 0.0822522*X.^2 + 0.480687*X + 0.0042848;
Y3_ref = 0.000128154*X.^9 - 0.00292649*X.^8 + 0.0281018*X.^7 - 0.14719*X.^6 + 0.455357*X.^5 - 0.839217*X.^4 + 0.884458*X.^3 - 0.529095*X.^2 + 0.368697*X + 0.00618416;

% Plot reference data with requested markers
plot(X, Y1_ref, 'k:', 'LineWidth', 2);
plot(X, Y2_ref, 'ks', 'MarkerSize', 6, 'MarkerIndices', round(linspace(1, length(X), 15)));
plot(X, Y3_ref, 'kd', 'MarkerSize', 6, 'MarkerIndices', round(linspace(1, length(X), 15)));

xlabel('time (hours)'); ylabel('Dimensionless Temperature (T*)');
% title('Stored Water (Fig 4)'); 
set(gca, 'Box', 'off', 'TickDir', 'in', 'LineWidth', 1.5, 'FontSize', 14);
legend({'Sim x/H=0.05', 'Sim x/H=0.55', 'Sim x/H=0.95', 'Ref x/H=0.05', 'Ref x/H=0.55', 'Ref x/H=0.95'}, 'Location', 'best', 'Box', 'off');
ylim([0 1]);
% ==============================================================
% Calculate and Print Correlation (R) and R-squared (R^2)
% ==============================================================
% For x/H = 0.05
R1_matrix = corrcoef(T_star(:, node_top), Y1_ref);
R1 = R1_matrix(1, 2);
R1_sq = R1^2;

% For x/H = 0.55
R2_matrix = corrcoef(T_star(:, node_mid), Y2_ref);
R2 = R2_matrix(1, 2);
R2_sq = R2^2;

% For x/H = 0.95
R3_matrix = corrcoef(T_star(:, node_bot), Y3_ref);
R3 = R3_matrix(1, 2);
R3_sq = R3^2;

fprintf('\n--- Correlation Analysis (Fig 4) ---\n');
MSE1 = mean((T_star(:, node_top) - Y1_ref).^2);
MSE2 = mean((T_star(:, node_mid) - Y2_ref).^2);
MSE3 = mean((T_star(:, node_bot) - Y3_ref).^2);

fprintf('x/H = 0.05 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R1, R1_sq, MSE1);
fprintf('x/H = 0.55 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R2, R2_sq, MSE2);
fprintf('x/H = 0.95 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R3, R3_sq, MSE3);
fprintf('------------------------------------\n');

% Plot Th* (Hot Coil Water - Fig 5)
subplot(1, 3, 2);
% Simulation plots
plot(t_hours, Th_star(:, node_top), 'k-', 'LineWidth', 1.5); hold on;
plot(t_hours, Th_star(:, node_mid), 'k--', 'LineWidth', 1.5);
plot(t_hours, Th_star(:, node_bot), 'k-.', 'LineWidth', 1.5);

% Reference polynomial plots for Fig 5
Y1_ref_fig5 = 0*X.^7 + 0.00167581*X.^6 - 0.0166634*X.^5 + 0.0598353*X.^4 - 0.0880729*X.^3 + 0.0254964*X.^2 + 0.0602337*X + 0.934334;
Y2_ref_fig5 = 0*X.^7 - 0.0030079*X.^6 + 0.0273408*X.^5 - 0.0968011*X.^4 + 0.184248*X.^3 - 0.268799*X.^2 + 0.396024*X + 0.500869;
Y3_ref_fig5 = 0*X.^7 - 0.000540285*X.^6 + 0.0013827*X.^5 + 0.0105175*X.^4 - 0.0368151*X.^3 - 0.0343903*X.^2 + 0.295194*X + 0.330674;

% Plot reference data with requested markers
plot(X, Y1_ref_fig5, 'k:', 'LineWidth', 2);
plot(X, Y2_ref_fig5, 'ks', 'MarkerSize', 6, 'MarkerIndices', round(linspace(1, length(X), 15)));
plot(X, Y3_ref_fig5, 'kd', 'MarkerSize', 6, 'MarkerIndices', round(linspace(1, length(X), 15)));

xlabel('time (hours)'); ylabel('Dimensionless Temperature (Th*)');
% title('Hot Coil Water (Fig 5)'); 
set(gca, 'Box', 'off', 'TickDir', 'in', 'LineWidth', 1.5, 'FontSize', 14);
legend({'Sim x/H=0.05', 'Sim x/H=0.55', 'Sim x/H=0.95', 'Ref x/H=0.05', 'Ref x/H=0.55', 'Ref x/H=0.95'}, 'Location', 'best', 'Box', 'off');
ylim([0 1]);
% ==============================================================
% Calculate and Print Correlation (R) and R-squared (R^2) for Fig 5
% ==============================================================
% For x/H = 0.05
R1_matrix_f5 = corrcoef(Th_star(:, node_top), Y1_ref_fig5);
R1_f5 = R1_matrix_f5(1, 2);
R1_sq_f5 = R1_f5^2;

% For x/H = 0.55
R2_matrix_f5 = corrcoef(Th_star(:, node_mid), Y2_ref_fig5);
R2_f5 = R2_matrix_f5(1, 2);
R2_sq_f5 = R2_f5^2;

% For x/H = 0.95
R3_matrix_f5 = corrcoef(Th_star(:, node_bot), Y3_ref_fig5);
R3_f5 = R3_matrix_f5(1, 2);
R3_sq_f5 = R3_f5^2;

fprintf('\n--- Correlation Analysis (Fig 5) ---\n');
MSE1_f5 = mean((Th_star(:, node_top) - Y1_ref_fig5).^2);
MSE2_f5 = mean((Th_star(:, node_mid) - Y2_ref_fig5).^2);
MSE3_f5 = mean((Th_star(:, node_bot) - Y3_ref_fig5).^2);

fprintf('x/H = 0.05 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R1_f5, R1_sq_f5, MSE1_f5);
fprintf('x/H = 0.55 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R2_f5, R2_sq_f5, MSE2_f5);
fprintf('x/H = 0.95 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R3_f5, R3_sq_f5, MSE3_f5);
fprintf('------------------------------------\n');

% Plot Tc* (Cold Coil Water - Fig 6)
subplot(1, 3, 3);
% Simulation plots
plot(t_hours, Tc_star(:, node_top), 'k-', 'LineWidth', 1.5); hold on;
plot(t_hours, Tc_star(:, node_mid), 'k--', 'LineWidth', 1.5);
plot(t_hours, Tc_star(:, node_bot), 'k-.', 'LineWidth', 1.5);

% Reference polynomial plots for Fig 6
Y1_ref_fig6 = 0*X.^7 - 0.00569982*X.^6 + 0.0566961*X.^5 - 0.225162*X.^4 + 0.479833*X.^3 - 0.686506*X.^2 + 0.83416*X - 0.0662528;
Y2_ref_fig6 = -0.000306233*X.^7 + 0.00454898*X.^6 - 0.0265726*X.^5 + 0.0810384*X.^4 - 0.137448*X.^3 + 0.0830661*X.^2 + 0.188162*X - 0.0109149;
Y3_ref_fig6 = 0.00343901*X.^7 - 0.0418129*X.^6 + 0.207598*X.^5 - 0.540516*X.^4 + 0.787415*X.^3 - 0.635301*X.^2 + 0.277486*X - 0.0404307;

% Plot reference data with requested markers
plot(X, Y1_ref_fig6, 'k:', 'LineWidth', 2);
plot(X, Y2_ref_fig6, 'ks', 'MarkerSize', 6, 'MarkerIndices', round(linspace(1, length(X), 15)));
plot(X, Y3_ref_fig6, 'kd', 'MarkerSize', 6, 'MarkerIndices', round(linspace(1, length(X), 15)));

xlabel('time (hours)'); ylabel('Dimensionless Temp (Tc*)');
% title('Cold Coil Water (Fig 6)'); 
set(gca, 'Box', 'off', 'TickDir', 'in', 'LineWidth', 1.5, 'FontSize', 14);
legend({'Sim x/H=0.05', 'Sim x/H=0.55', 'Sim x/H=0.95', 'Ref x/H=0.05', 'Ref x/H=0.55', 'Ref x/H=0.95'}, 'Location', 'best', 'Box', 'off');
ylim([0 1]);
% ==============================================================
% Calculate and Print Correlation (R), R-squared (R^2), and MSE for Fig 6
% ==============================================================
% For x/H = 0.05
R1_matrix_f6 = corrcoef(Tc_star(:, node_top), Y1_ref_fig6);
R1_f6 = R1_matrix_f6(1, 2);
R1_sq_f6 = R1_f6^2;
MSE1_f6 = mean((Tc_star(:, node_top) - Y1_ref_fig6).^2);

% For x/H = 0.55
R2_matrix_f6 = corrcoef(Tc_star(:, node_mid), Y2_ref_fig6);
R2_f6 = R2_matrix_f6(1, 2);
R2_sq_f6 = R2_f6^2;
MSE2_f6 = mean((Tc_star(:, node_mid) - Y2_ref_fig6).^2);

% For x/H = 0.95
R3_matrix_f6 = corrcoef(Tc_star(:, node_bot), Y3_ref_fig6);
R3_f6 = R3_matrix_f6(1, 2);
R3_sq_f6 = R3_f6^2;
MSE3_f6 = mean((Tc_star(:, node_bot) - Y3_ref_fig6).^2);

fprintf('\n--- Correlation Analysis (Fig 6) ---\n');
fprintf('x/H = 0.05 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R1_f6, R1_sq_f6, MSE1_f6);
fprintf('x/H = 0.55 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R2_f6, R2_sq_f6, MSE2_f6);
fprintf('x/H = 0.95 : Correlation (R) = %f, R^2 = %f, MSE = %e\n', R3_f6, R3_sq_f6, MSE3_f6);
fprintf('------------------------------------\n');