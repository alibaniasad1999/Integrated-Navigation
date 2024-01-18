%% model aided filtering
clear;
clc;
%% load data
load('Eval.mat');
%%
ax = a(:,1);
ay = a(:,2);
az = a(:,3);

theta = atan2(ax,sqrt(ay.^2+az.^2));
theta = [time', theta];
u = [time', u];



