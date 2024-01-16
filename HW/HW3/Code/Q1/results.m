clear;
clc;
load('DataSet_02.mat')
dt = 0.01;
position_num = in_profile(1, 2:4);
velocity_num = in_profile(1, 5:7);
attitude_num = in_profile(1, 8:10);
%% initial state
x_num = [position_num, velocity_num, attitude_num]';
u_num = IMU_true(1, 2:7)';
position_array = zeros(length(IMU_true), 3);
velocity_array = zeros(length(IMU_true), 3);
attitude_array = zeros(length(IMU_true), 3);
out_profile = zeros(length(IMU_true), 9);
u_old = 0;
x_old = 0;
delta_x_num = zeros(9, 1);
%% run simulation
for i = 1:length(IMU_true)
    u_num = IMU_true(i, 2:7)';
    A_num = A_func(x_num, u_num);
    B_num = B_func(x_num, u_num);
    delta_u_num = u_num - u_old;
    delta_x_dot_num = A_num*delta_x_num + B_num*delta_u_num;
    delta_x_num = delta_x_num + delta_x_dot_num*dt;
    x_num = x_num + delta_x_num;
    u_old = u_num;
    out_profile(i, :) = x_num';
end
%% plot

figure(1)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 2), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 1), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('latitude', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/latitude_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(2)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 3), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 2), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('longitude', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/longitude_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(3)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 4), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 3), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('altitude', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/altitude_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(4)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 5), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 4), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$v_x$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/velocity_x_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(5)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 6), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 5), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$v_y$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/velocity_y_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(6)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 7), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 6), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$v_z$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/velocity_z_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(7)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 8), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 7), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\phi$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/phi_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(8)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 9), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 8), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\theta$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/theta_lin','-depsc');
else
    fprintf("Ooooooops\n")
end

figure(9)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 10), 'LineWidth', 2, 'Color','r');
plot(IMU_true(:, 1), out_profile(:, 9), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\psi$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight
[dir_state, ~, ~] = mkdir('../../../Figure/Q5');
if dir_state
    print('../../../Figure/Q5/psi_lin','-depsc');
else
    fprintf("Ooooooops\n")
end
