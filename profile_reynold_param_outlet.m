% Profile for reynold_param_outlet.m (Transient)
addpath('..');
model = model_data();
tank = tank_data();
N = model.N;
Reh_array = [2.17e3, 1.73e4, 4.32e4, 6.48e4]; % Rec values in this script
t_span = linspace(0, 3*3600, 200);

figure('Position', [100, 100, 1200, 400]);

for k = 1:length(Reh_array)
    disp(['Running simulation for Rec = ', num2str(Reh_array(k))]);
    model.Rec = Reh_array(k);
    model.Reh = 1.62e5; % From original script
    T_in = ones(N, 1) * model.T0;
    Th_in = ones(N, 1) * model.T0;
    Tc_in = ones(N, 1) * model.T0;
    T_initial = [T_in; Th_in; Tc_in];
    
    [t_out, T_out] = ode15s(@(t, T_state) heat_trans_model(t, T_state, model, tank), t_span, T_initial);
    T_raw = T_out(:, 1:N);
    T_star = (T_raw - 300) ./ 100;
    
    t_hours = t_out / 3600;
    x_H = (0.5 : 1 : N) / N;
    
    subplot(1, length(Reh_array), k);
    [X, Y] = meshgrid(t_hours, x_H);
    contourf(X, Y, T_star', 100, 'LineColor', 'none');
    set(gca, 'YDir', 'reverse');
    colormap(jet);
    caxis([0 1]);
    
    xlabel('Time (Hours)', 'FontSize', 12);
    if k == 1
        ylabel('Depth (x/H)', 'FontSize', 12);
    end
    title(sprintf('Re_c = %1.2e', Reh_array(k)), 'FontSize', 12);
    xticks(0:0.5:3);
end

% Add a shared colorbar
h = axes('Position', [0 0 1 1], 'Visible', 'off');
c = colorbar(h, 'Position', [0.93 0.15 0.02 0.75]);
caxis(h, [0 1]);
colormap(h, jet);
ylabel(c, 'Dimensionless Temperature (T^*)', 'FontSize', 12);

saveas(gcf, 'profile_reynold_param_outlet.png');
disp('Saved profile_reynold_param_outlet.png');
