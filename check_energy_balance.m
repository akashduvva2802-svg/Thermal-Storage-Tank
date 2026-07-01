model = model_data();
tank = tank_data();
model.Reh_list = [2.17e3, 6.48e4, 1.62e5, 2.43e5];
Reh_array = model.Reh_list;
t_span = linspace(0, 3600*10, 200);

for k=1:length(Reh_array)
    model.Reh = Reh_array(k);
    model.Rec = 1.73e4;
    T_in = ones(model.N, 1)*300;
    T_initial = [T_in; T_in; T_in];
    [t_out, T_out] = ode15s(@(t, T_state) heat_trans_model(t, T_state, model, tank), t_span, T_initial);
    T_ss = T_out(end, :)';
    
    T = T_ss(1:model.N);
    Th = T_ss(model.N+1:2*model.N);
    Tc = T_ss(2*model.N+1:3*model.N);
    
    prop_in_h = water_properties(model.Th_in);
    prop_in_c = water_properties(model.Tc_in);
    mh = model.Reh*(pi*tank.di*prop_in_h.mu)/4;
    mc = model.Rec*(pi*tank.di*prop_in_c.mu)/4;
    
    % Sum over all nodes for Q_hot and Q_cold as implemented in calc_Q_hot
    delta = tank.di / tank.Dcoil;
    Ao_h = pi*tank.do*(tank.Lcoil_h * model.dx / tank.H);
    Ao_c = pi*tank.do*(tank.Lcoil_c * model.dx / tank.H);
    
    Q_hot_total = 0;
    Q_cold_total = 0;
    
    for i = 1:model.N
        % Hot
        Tsh_guess = (Th(i) + T(i))/2;
        for iter = 1:3
            [UA_h, ho_h] = coil_overall_htc(Th(i), T(i), Tsh_guess, delta, 'hot', model);
            if ho_h == 0
                Tsh_guess = Th(i);
            else
                Tsh_guess = T(i) + (UA_h/(ho_h*Ao_h))*(Th(i) - T(i));
            end
        end
        Q_hot_total = Q_hot_total + UA_h*(Th(i) - T(i));
        
        % Cold
        Tsc_guess = (Tc(i) + T(i))/2;
        for iter = 1:3
            [UA_c, ho_c] = coil_overall_htc(Tc(i), T(i), Tsc_guess, delta, 'cold', model);
            if ho_c == 0
                Tsc_guess = Tc(i);
            else
                Tsc_guess = T(i) + (UA_c/(ho_c*Ao_c))*(Tc(i) - T(i));
            end
        end
        Q_cold_total = Q_cold_total + UA_c*(T(i) - Tc(i));
    end
    
    Q_hot_in_fluid = mh*prop_in_h.Cp*(model.Th_in - Th(end));
    Q_cold_out_fluid = mc*prop_in_c.Cp*(Tc(1) - model.Tc_in);
    
    Q_loss = 0;
    for i=1:model.N
        x=(model.N-i+0.5)*model.dx;
        Ts_guess=(T(i)+300)/2;
        for iter=1:3
            [UA_loss, hamb]=loss_htc(Ts_guess,300,x);
            if hamb~=0
                Ts_guess=300+(UA_loss/(hamb*pi*tank.D*model.dx))*(T(i)-300);
            end
        end
        Q_loss = Q_loss + UA_loss*(T(i)-300);
    end
    
    fprintf('Re_h=%.2e\n', Reh_array(k));
    fprintf('  Fluid enthalpy drop  : Q_hot=%.2f, Q_cold=%.2f\n', Q_hot_in_fluid, Q_cold_out_fluid);
    fprintf('  Tank heat transfer   : Q_hot=%.2f, Q_cold=%.2f\n', Q_hot_total, Q_cold_total);
    fprintf('  Heat Loss            : %.2f\n', Q_loss);
    fprintf('  Balance (Tank Q)     : Q_hot - Q_cold - Q_loss = %.2f\n\n', Q_hot_total - Q_cold_total - Q_loss);
end
