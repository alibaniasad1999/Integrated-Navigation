clc;
clear;
%% run section a %%
a;
%% checking observability %%
% check observability of the system
if rank(obsv(A,C)) == size(A,1)
    disp('System is observable');
else
    disp('System is not observable');
end