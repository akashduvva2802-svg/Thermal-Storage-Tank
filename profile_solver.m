% Profile for solver.m (Transient)
addpath('..');
t_span = linspace(0, 3*3600, 200); % 3 hours
model = model_data();
tank = tank_data();
N = model.N;
T_in = ones(N, 1) * model.T0;
Th_in = ones(N, 1) * model.T0;
Tc_in = ones(N, 1) * model.T0;
T_initial = [T_in; Th_in; Tc_in];

[t_out, T_out] = ode15s(@(t, T_state) heat_trans_model(t, T_state, model, tank), t_span, T_initial);

T_raw  = T_out(:, 1:N);
T_t0 = 300; delta_T_in = 100;
T_star = (T_raw - T_t0) ./ delta_T_in;

t_hours = t_out / 3600;
x_H = (0.5 : 1 : N) / N;

figure('Position', [100, 100, 600, 500]);
[X, Y] = meshgrid(t_hours, x_H);
contourf(X, Y, T_star', 100, 'LineColor', 'none');
set(gca, 'YDir', 'reverse');
colormap(jet);
caxis([0 1]);
h_cb = colorbar;
ylabel(h_cb, 'Dimensionless Temperature (T^*)');
xlabel('Time (Hours)', 'FontSize', 12);
ylabel('Depth (x/H)', 'FontSize', 12);
title('Transient Profile (solver)', 'FontSize', 14);
xticks(0:0.25:3);

saveas(gcf, 'profile_solver.png');
disp('Saved profile_solver.png');
