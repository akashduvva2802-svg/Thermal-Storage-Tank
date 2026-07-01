% Profile for steady_param_hot.m (Steady State)
addpath('..');
model = model_data();
tank = tank_data();
N = model.N;
Reh_array = [2.16e3, 6.48e4, 1.6e5, 2.4e5];
t_span = linspace(0, 3*3600, 200);

figure('Position', [100, 100, 1000, 400]);

for k = 1:length(Reh_array)
    disp(['Running simulation for Reh = ', num2str(Reh_array(k))]);
    model.Reh = Reh_array(k);
    T_in = ones(N, 1) * model.T0;
    Th_in = ones(N, 1) * model.T0;
    Tc_in = ones(N, 1) * model.T0;
    T_initial = [T_in; Th_in; Tc_in];
    
    [~, T_out] = ode15s(@(t, T_state) heat_trans_model(t, T_state, model, tank), t_span, T_initial);
    T_steady = T_out(end, 1:N);
    T_star_steady = (T_steady - 300) ./ 100;
    
    subplot(1, length(Reh_array), k);
    x_H = (0.5 : 1 : N) / N;
    T_grid = repmat(T_star_steady', 1, 2);
    
    imagesc([0 1], x_H, T_grid);
    set(gca, 'YDir', 'reverse');
    colormap(jet);
    caxis([0 1]);
    
    xlabel('Steady State', 'FontSize', 12);
    if k == 1
        ylabel('Depth (x/H)', 'FontSize', 12);
    end
    title(sprintf('Re_h = %1.2e', Reh_array(k)), 'FontSize', 12);
    set(gca, 'XTick', []);
end

% Add a shared colorbar
h = axes('Position', [0 0 1 1], 'Visible', 'off');
c = colorbar(h, 'Position', [0.93 0.15 0.02 0.75]);
caxis(h, [0 1]);
colormap(h, jet);
ylabel(c, 'Dimensionless Temperature (T^*)', 'FontSize', 12);

saveas(gcf, 'profile_steady_param_hot.png');
disp('Saved profile_steady_param_hot.png');
