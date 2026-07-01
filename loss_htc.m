function [UA_loss, h_amb] = loss_htc(Ts, Tinf, x)
model = model_data();
tank = tank_data();
Ro = tank.D/2;
dx = model.dx;
Ao = tank.D*pi*dx;
Nu_amb = ambient_nu(Ts, Tinf, x);
h_amb = Nu_amb*tank.kair/tank.D;
    % Override UA_loss with the constant value from the paper
    UA_loss = 0.2026 / model.N;
end

function Gr = grashoff(Ts, Tinf, x)
tank = tank_data();
g = 9.81;
Gr = g*(Ts - Tinf)*(x^3)/((tank.nu_air^2)*(Ts + Tinf));
end

function Nu_amb = ambient_nu(Ts, Tinf, x)
model = model_data();
tank = tank_data();
Gr = grashoff(Ts, Tinf, x);
Pr = tank.Pr_air;
Nu_amb = ((7*Gr*Pr^2)/(100 + 105*Pr))^0.25 + 4*(272 + Pr)*x/(35*(64 + 63*Pr)*tank.D);
end
