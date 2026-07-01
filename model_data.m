function model = model_data()
tank = tank_data();
model.Th_in = 400;
model.Tc_in = 300;
model.T0 = 300;
model.Tinf = 300;
model.Reh = 1.62e5;
model.Rec = 1.73e04;
model.N = 10;
model.H = tank.H;
model.dx = model.H/model.N;
model.dt = 10;
end
