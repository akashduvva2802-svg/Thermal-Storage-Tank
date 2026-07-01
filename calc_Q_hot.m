function Q_hot_kW = calc_Q_hot(T_out, model, tank)
    tank = tank_data();
    N = model.N;
    Q_hot_kW = zeros(size(T_out, 1), 1);
    
    Tinf = model.T0;
    delta = tank.di / tank.Dcoil;
    Ao = pi*tank.do*(tank.Lcoil_h * model.dx / tank.H);
    
    for j = 1:size(T_out, 1)
        T_state = T_out(j, :)';
        T = T_state(1:N);
        Th = T_state(N+1:2*N);
        
        Q_total = 0;
        
        for i = 1:N
            % Fixed-Point Iteration for Surface Temperatures
            Tsh_guess = (Th(i) + T(i))/2;
            for iter = 1:3
                [UA_h, ho_h] = coil_overall_htc(Th(i), T(i), Tsh_guess, delta, 'hot', model);
                if ho_h == 0
                    Tsh_guess = Th(i);
                else
                    Tsh_guess = T(i) + (UA_h/(ho_h*Ao))*(Th(i) - T(i));
                end
            end
            
            Q_node = UA_h * (Th(i) - T(i));
            Q_total = Q_total + Q_node;
        end
        Q_hot_kW(j) = Q_total / 1000; % Convert W to kW
    end
end
