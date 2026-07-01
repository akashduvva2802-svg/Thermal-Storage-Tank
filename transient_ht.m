function dTdt_state = transient_ht(t, T_state, model, tank)
omega = 1/1800; % as taken in paper
N = model.N;

T = T_state(1:N);
Th = T_state(N+1:2*N);
Tc = T_state(2*N + 1:3*N);

dTdt = zeros(N,1);
dThdt = zeros(N,1);
dTcdt = zeros(N,1);

% Initialize surface temperatures
Ts = zeros(N, 1);
Tsh = zeros(N, 1);
Tsc = zeros(N, 1);

gamma = 1; % for now I consider that coil is throughout

Tinf = model.T0; 

Ac = pi*(tank.D^2)/4;
A = pi*tank.D*model.dx; % Outer surface area of a tank node
delta = tank.di / tank.Dcoil;
Ao = pi*tank.do*(tank.Lcoil_h * model.dx / tank.H); % Outer surface area of coil in one node

prop_in_h = water_properties(function_Th_in(t, omega));
prop_in_c = water_properties(model.Tc_in);
mh = model.Reh*(pi*tank.di*prop_in_h.mu)/4;
mc = model.Rec*(pi*tank.di*prop_in_c.mu)/4;

for i = 1:N
   if i == 1 
       Th_enter = function_Th_in(t,omega);
   else
       Th_enter = Th(i-1); 
   end
   if i == N
       Tc_enter = model.Tc_in;
   else
       Tc_enter = Tc(i+1); 
   end
   
   if i == 1 
       T_above = T(i);  
   else
       T_above = T(i-1); 
   end
   
   if i == N
       T_below = T(i);  
   else
       T_below = T(i+1); 
   end
   
   prop = water_properties(T(i));
   prop_h = water_properties(Th(i));
   prop_c = water_properties(Tc(i));
   
   % Changed do to di for the mass flow calculations
   m = (pi*(tank.D^2/4))*model.dx*prop.rho;
   s = tank.Lcoil_h * model.dx / tank.H;
   m_node_h = prop_h.rho * s * pi*(tank.di^2)/4;
   m_node_c = prop_c.rho * s * pi*(tank.di^2)/4;
   % --- INNER LOOP: Fixed-Point Iteration for Surface Temperatures ---
   % Initial Guess of Ts
   Tsh_guess = (Th(i) + T(i))/2;
   Tsc_guess = (Tc(i) + T(i))/2;
   Ts_guess = (T(i) + Tinf)/2;
   
   x = (N-i + 0.5)*model.dx;
   
   % Iterate to convergence
   for iter = 1:5
       [UA_h, ho_h] = coil_overall_htc(Th(i), T(i), Tsh_guess, delta, 'hot', model);
       [UA_c, ho_c] = coil_overall_htc(Tc(i), T(i), Tsc_guess, delta, 'cold', model); 
       [UA_loss, h_amb] = loss_htc(Ts_guess, Tinf, x);

       % Update the guesses using the newly calculated HTCs
       if ho_h == 0
           Tsh_guess = Th(i);
       else
           Tsh_guess = T(i) + (UA_h/(ho_h*Ao))*(Th(i) - T(i));
       end
       
       if ho_c == 0
           Tsc_guess = Tc(i);
       else
           Tsc_guess = T(i) + (UA_c/(ho_c*Ao))*(Tc(i) - T(i));
       end
       
       if h_amb == 0
           Ts_guess = T(i);
       else
           Ts_guess = Tinf + (UA_loss/(h_amb*A))*(T(i) - Tinf);
       end
   end
   
   % converged values of Ts
   Tsh(i) = Tsh_guess;
   Tsc(i) = Tsc_guess;
   Ts(i) = Ts_guess;
   % ------------------------------------------------------------------

   dTdt(i) = (1/(m*prop.Cp))*(gamma*UA_h*(Th(i) - T(i)) ...
             + prop.k*Ac*(T_above - T(i))/model.dx  ...
             + prop.k*Ac*(T_below - T(i))/model.dx ...
             + gamma*UA_c*(Tc(i) - T(i)) ...
             - UA_loss*(T(i) - Tinf));
             
   dThdt(i) = (1/(m_node_h*prop_h.Cp))*(mh*prop_h.Cp*(Th_enter - Th(i)) - gamma*UA_h*(Th(i) - T(i)));
   dTcdt(i) = (1/(m_node_c*prop_c.Cp))*(- mc*prop_c.Cp*(Tc(i) - Tc_enter) + gamma*UA_c*(T(i) - Tc(i)));
end

dTdt_state = [dTdt; dThdt; dTcdt];

end
