%% init 
load('Q2.mat')
theta_m = atan2(-ax, sqrt(ay.^2 + az.^2));
theta_sim = [time', theta_m];
u_sim = [time', u];
g = [time', g];