function prop = properties(T)
prop.mu = 79.86*exp(-0.04086*T) + 0.005841*exp(-0.008327*T);
prop.k = (-8.308e-06)*(T^2) + 0.0065*T - 0.6098;
prop.Cp = 0.01259*(T^2) - 8.035*T + 5460;
prop.rho = -0.002934*(T^2) + 1.47*T + 819.2;
prop.nu = prop.mu/prop.rho;
end
