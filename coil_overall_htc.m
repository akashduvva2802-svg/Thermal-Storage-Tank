function [UA, ho] = coil_overall_htc(T, Tinf, Ts, delta, type, model)
    % type: 'hot' or 'cold'
    if nargin < 5
        type = 'hot';
    end
    tank = tank_data();
    % Fluid properties inside the tube
    prop = water_properties(T);
    % Fluid properties at the bulk temperature (for external convection)
    prop1 = water_properties(Tinf); % Evaluate at bulk temperature
    
    [Nu_d, Nu_l] = nusselt(T, delta, Ts, Tinf, type, model);
    hi = Nu_d*prop.k/tank.di;
    
    if strcmp(type, 'hot')
        Lcoil = tank.Lcoil_h;
    else
        Lcoil = tank.Lcoil_c;
    end
    s = Lcoil * (model.dx / tank.H); % actual length of coil in one node
    Ai = pi*tank.di*s;
    Ao = pi*tank.do*s;

    ho = Nu_l*prop1.k/tank.do;
    % Calculate overall heat transfer coefficient UA
    UA = 1/((1/(hi*Ai)) + (1/(ho*Ao)) + (log(tank.do/tank.di)/(2*pi*tank.khx*s)));
end

function f = frict_fac(Re, delta)
    if delta >= 0.034 && delta <= 300
        f = 0.304*Re.^(-0.25) + 0.029*(delta^0.5);
    else 
        disp('error_delta')
        f = NaN;
    end
end

function Re = reynold(T, type, model)
    prop = water_properties(T);
    if strcmp(type, 'hot')
        Re0 = model.Reh;
        prop_in = water_properties(model.Th_in);
    else
        Re0 = model.Rec;
        prop_in = water_properties(model.Tc_in);
    end
    Re = Re0;
end

function Gr = grashoff(Ts, T_inf, x, nu)
    g = 9.81;
    prop_s = water_properties(Ts);
    prop_inf = water_properties(T_inf);
    beta_dT = abs(prop_inf.rho - prop_s.rho) / prop_inf.rho;
    Gr = (g * beta_dT / (nu^2)) * x^3;
end

function [Nu_D, Nu_L] = nusselt(T, delta, Ts, T_inf, type, model)
    prop = water_properties(T);
    Pr = (prop.mu)*(prop.Cp)/(prop.k);
    
    if nargin < 5
        type = 'hot';
    end
    
    Re = reynold(T, type, model);
    
    % Eq (13) Critical Reynolds number
    Re_cr = 2.1e3 * (1 + delta);
    
    % Inner Nusselt number Nu_D
    if Re > Re_cr
        % Turbulent flow - Pethukov correlation Eq (10)
        f = frict_fac(Re, delta);
        Nu_D = (f/8) * Re.* Pr / (1.07 + 12.7 * (f/8)^0.5 * (Pr^(0.667) - 1));
    else
        % Laminar flow - Xin and Ebadian correlation Eq (12)
        % (Corrected exponent from paper typo: Re^2 -> Re^0.92)
        if strcmp(type, 'hot')
            % For cooling of the fluid (hot coil), n=0.3
            Nu_D = 0.00619 * Re.^(0.92) * (Pr^0.3) * (1 + 3.455 * delta);
        else
            % For heating of the fluid (cold coil), n=0.4
            Nu_D = 0.00619 * Re.^(0.92) * (Pr^0.4) * (1 + 3.455 * delta);
        end
    end
    
    % Outer Nusselt number Nu_L - Ali correlation Eq (14)
    if nargin >= 4
        % Hardcode characteristic length to outer tube diameter
        tank = tank_data();
        L_char = tank.do; 
        
        prop_f = water_properties(T_inf); % properties evaluated at bulk temp instead of film
        nu_f = prop_f.nu;
        Pr_f = prop_f.mu * prop_f.Cp / prop_f.k;
        
        Gr_L = grashoff(Ts, T_inf, L_char, nu_f);
        Ra_L = Gr_L * Pr_f;
        Nu_L = 0.106 * abs(Ra_L)^0.335;
    else
        Nu_L = NaN;
    end
end
