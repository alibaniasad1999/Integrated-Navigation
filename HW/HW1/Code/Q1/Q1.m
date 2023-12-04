%% HW1 Q1 %%
%% a ---> IV %%
load('DataSet_01.mat')
set(gca, 'FontSize', 16)
plot(IMU_true(:, 1), IMU_true(:, 2:4), 'linewidth', 2);
legend('a_x', 'a_x', 'a_x', 'Location','southeast', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('acceleration($m/s^2$)', 'interpreter', 'latex', 'FontSize', 24);
title('');
[dir_state, ~, ~] = mkdir('../../Figure/Q1');
if dir_state
    print('../../Figure/Q1/immovability','-depsc');
else
    fprintf("Ooooooops\n")
end
%% a ----> V %%
% scenario 1
% stationary
load('DataSet_01.mat')
roll = atan2(IMU_true(:, 3), IMU_true(:, 4));
pitch = atan2(-IMU_true(:, 2), sqrt(IMU_true(:, 3).^2 + IMU_true(:, 4).^2));

set(gca, 'FontSize', 16)
plot(IMU_true(:, 1), [roll, pitch], 'linewidth', 2);
legend('roll', 'pitch', 'Location','southeast', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('angle($rad$)', 'interpreter', 'latex', 'FontSize', 24);
title('');
[dir_state, ~, ~] = mkdir('../../Figure/Q1');
if dir_state
    print('../../Figure/Q1/roll_pitch','-depsc');
else
    fprintf("Ooooooops\n")
end
% scenario 2 
% moving 
load('DataSet_02.mat')
roll = atan2(IMU_true(:, 3), IMU_true(:, 4));
pitch = atan2(-IMU_true(:, 2), sqrt(IMU_true(:, 3).^2 + IMU_true(:, 4).^2));

set(gca, 'FontSize', 16)
plot(IMU_true(:, 1), [roll, pitch], 'linewidth', 2);
legend('roll', 'pitch', 'Location','southeast', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('angle($rad$)', 'interpreter', 'latex', 'FontSize', 24);
title('');
[dir_state, ~, ~] = mkdir('../../Figure/Q1');
if dir_state
    print('../../Figure/Q1/roll_pitch_2','-depsc');
else
    fprintf("Ooooooops\n")
end

%% b ----> IV %%
load('DataSet_01.mat')
set(gca, 'FontSize', 16)
plot(IMU_true(:, 1), IMU_true(:, 5:7), 'linewidth', 2);
legend('w_x', 'w_x', 'w_x', 'Location','southeast', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('angular velocity($rad/s$)', 'interpreter', 'latex', 'FontSize', 24);

[dir_state, ~, ~] = mkdir('../../Figure/Q1');
if dir_state
    print('../../Figure/Q1/immovability_w','-depsc');
else
    fprintf("Ooooooops\n")
end

%% c ----> b %%
% scenario 1
load('DataSet_01.mat')
m = MAG_true(:, 2:4);
a = IMU_true(:, 2:4);
p = cross(m, a);
psi = (p(:, 1)./(a(:, 2).*p(:, 3) - a(:, 3).*p(:, 2)));

set(gca, 'FontSize', 16)
plot(MAG_true(:, 1), psi, 'linewidth', 2);
legend('psi', 'Location','northeast', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('angle($rad$)', 'interpreter', 'latex', 'FontSize', 24);
title('');
[dir_state, ~, ~] = mkdir('../../Figure/Q1');
if dir_state
    print('../../Figure/Q1/psi','-depsc');
else
    fprintf("Ooooooops\n")
end

% scenario 2
load('DataSet_02.mat')
m = MAG_true(:, 2:4);
a = IMU_true(:, 2:4);
p = cross(m, a);
psi = (p(:, 1)./(a(:, 2).*p(:, 3) - a(:, 3).*p(:, 2)));

set(gca, 'FontSize', 16)
plot(MAG_true(:, 1), psi, 'linewidth', 2);
legend('psi', 'Location','northeast', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('angle($rad$)', 'interpreter', 'latex', 'FontSize', 24);
title('');
[dir_state, ~, ~] = mkdir('../../Figure/Q1');
if dir_state
    print('../../Figure/Q1/psi_2','-depsc');
else
    fprintf("Ooooooops\n")
end



