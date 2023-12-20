clc;
clear;
%% linearization of INS equations using MATLAB symbolic toolbox
syms phi theta psi v_n v_e v_d L lambda h real
syms w_x w_y w_z a_x a_y a_z real
% constants
syms g R_m R_p R_g Omega_e real

% state vector
x = [phi theta psi v_n v_e v_d L lambda h]';
% input vector
u = [w_x w_y w_z a_x a_y a_z]';
% state vector derivative
omega_ie_n = [Omega_e*cos(L); 0; -Omega_e*sin(L)];
omega_en_n = [v_e/(R_p+h); -v_n/(R_m+h); -v_e*tan(L)/(R_p+h)];
omega_in_n = omega_ie_n + omega_en_n;
C_nb = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi),...
cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi),...
sin(phi)*cos(theta);
sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)...
, -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi),...
cos(phi)*cos(theta)];
C_bn = C_nb';
omega_nb_b = [w_x; w_y; w_z] - C_bn'*omega_in_n;
attitude_rate = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
    0 cos(phi) -sin(phi);
    0 sin(phi)/cos(theta) cos(phi)/cos(theta)]*omega_nb_b;
g_n = [0; 0; g];
g_n = g_n*(1 + h/R_g)^2;
F_n = C_bn * [a_x; a_y; a_z];
correction_term = cross(2*omega_ie_n + omega_en_n, [v_n; v_e; v_d]);
velocity_rate = F_n + g_n - correction_term;
position_rate = [v_n/(R_m+h); v_e/((R_p+h)*cos(L)); -v_d];
state_rate = [attitude_rate; velocity_rate; position_rate];
% linearization
A = jacobian(state_rate, x);
B = jacobian(state_rate, u);
%% save to file
matlabFunction(A, 'File', 'A_matrix', 'Vars', {x, u, g, R_m, R_p, R_g, Omega_e});
matlabFunction(B, 'File', 'B_matrix', 'Vars', {x, u, g, R_m, R_p, R_g, Omega_e});
%% test in zero
x = zeros(9, 1);
u = zeros(6, 1);
g = 9.7803267714;
R_m = 6378137;
R_p = 6356752.3142;
R_g = 6371008.7714;
Omega_e = 7.292116e-5;
A = A_matrix(x, u, g, R_m, R_p, R_g, Omega_e);
B = B_matrix(x, u, g, R_m, R_p, R_g, Omega_e);

%% load data
load('DataSet_02.mat');
%% initial conditions
position = in_profile(1, 2:4)';
velocity = in_profile(1, 5:7)';
attitude = in_profile(1, 8:10)';
%% parameters:
%% Earth Parameters
f = 1/298.2572236;
e = 0.0818191908426;
dt = 0.01;
R_0 = 6378137;
ee = e*e;
Omega_e = 7.29211651452621e-5; %Earth rotation rate
a = 6378178;
R_m = a*(1-ee)/((1-ee*(sin(position(1)))^2)^(3/2));
R_p = a/((1-ee*(sin(position(1)))^2)^(1/2));
R_G = sqrt(R_m*R_p);
%% initial state
x = [attitude; velocity; position];
%% initial input
u = [IMU_true(1, 5:7)'; IMU_true(1, 2:4)'];
%% use symbolic function
out_profile = zeros(length(IMU_true), 9);
for i = 1:length(IMU_true)
    u = [IMU_true(i, 5:7)'; IMU_true(i, 2:4)'];
    attitude = x(1:3);
    velocity = x(4:6);
    position = x(7:9);
    out_profile(i, :) = [position', velocity', attitude'];
    A = A_matrix(x, u, g, R_m, R_p, R_G, Omega_e);
    B = B_matrix(x, u, g, R_m, R_p, R_G, Omega_e);
    x = x + A*x*dt + B*u*dt;
end
%% plot
figure(1);
plot(IMU_true(:, 1), out_profile(:, 7), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 8), 'b');
hold off;

figure(2);
plot(IMU_true(:, 1), out_profile(:, 8), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 9), 'b');
hold off;

figure(3);
plot(IMU_true(:, 1), out_profile(:, 9), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 10), 'b');
hold off;

figure(4);
plot(IMU_true(:, 1), out_profile(:, 4), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 5), 'b');
hold off;

figure(5);
plot(IMU_true(:, 1), out_profile(:, 5), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 6), 'b');
hold off;

figure(6);
plot(IMU_true(:, 1), out_profile(:, 6), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 7), 'b');
hold off;

figure(7);
plot(IMU_true(:, 1), out_profile(:, 1), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 2), 'b');
hold off;

figure(8);
plot(IMU_true(:, 1), out_profile(:, 2), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 3), 'b');
hold off;

figure(9);
plot(IMU_true(:, 1), out_profile(:, 3), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 4), 'b');
hold off;



