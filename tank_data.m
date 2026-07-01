function tank = tank_data()
tank.khx = 30;
tank.p = 36.2e-03;
tank.Dcoil = 0.49;
tank.di = 21.6e-03;
tank.do = 26.9e-03;
tank.Lcoil_h = 85.1/2; % Length of hot coil in meters
tank.Lcoil_c = 85.1/2; % Length of cold coil in meters
tank.H = 2;
tank.D = 1.25;
tank.kmat = 0.0474; 
tank.dins = 200e-03; 
tank.kair = 0.026; %W/mK at 25°C
tank.Pr_air = 0.71; % at 25°C
tank.nu_air = 1.56e-05; % at 25°C
end
