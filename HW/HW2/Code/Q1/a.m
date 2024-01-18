clc;
clear;
%% state space %%
syms x1 x2 x3 x4 x5 u
% states 1 ---> height
% states 2 ---> velocity
% states 3 ---> bias
% states 4 ---> a0
% states 5 ---> a1
% constant
h_r = 0.013;
R = 0.0333;
m = 0.063;
%% state space
f1 = x2;
z = x1 + x3;
h_bar = (-z + h_r)/(2.5*R);
K_G = x4 + x5*(2/h_bar)^2;
f2 = -9.81 - K_G/m*u;
f3 = 0;
f4 = 0;
f5 = 0;
f = [f1;f2;f3;f4;f5];
%% jacobi matrix
A = jacobian(f,[x1,x2,x3,x4,x5]);
B = jacobian(f,u);
%% linearization
A = subs(A,[x1,x2,x3,x4,x5,u],[0,0,0.2,0.03794,0.9621,m*9.81]);
B = subs(B,[x1,x2,x3,x4,x5,u],[0,0,0.2,0.03794,0.9621,m*9.81]);
A = double(vpa(A));
% B = double(vpa(B));
C = [1,0,0,0,0];