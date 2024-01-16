clc
close all
%%
tf=20;
Ts=0.01;
x = [-2;0;0.2;0.9621;0.03794];
P = 1e-2*eye(5);
R=1;
Q=0.002;
xhat=x;
%%
N=100;
xpart = zeros(5, N);
% Initialize the particle filter.
for i = 1 : N
    xpart(:, i) = x + sqrt(P) * randn(5, 1);
end
%%
tnew=[0];
xnew=[x];
xhatnew=[xhat];
P11new=[P(1,1)];
xhatpfnew=[x];
xpartminus = zeros(5, N);
%%
out = sim('Quad_IGE_Landing');
%% load data u and z %%
%% run simulink %%
u_in = out.u.Data;
zm = out.zm.Data; % measurement
for t=1:length(zm)
% %% System 
%    x=0.5*x+25*x/(1+x^2)+8*cos(1.2*(t))+sqrt(Q)*randn;
%    ym=x^2/20+sqrt(R)*randn;
% %% kf
% xhat=0.5*xhat+25*xhat/(1+xhat^2)+8*cos(1.2*t);
% F=0.5+25*(1-xhat^2)/((1+xhat^2)^2);
% P=F*P*F'+Q;
% C=xhat/10;
% K=P*C'*inv(C*P*C'+R);
% yhat=xhat^2/20;
% xhat=xhat+K*(ym-yhat);
% P=(eye(1)-K*C)*P*(eye(1)-K*C)'+K*R*K';
%% PF
% Particle filter
for i = 1 : N
    x1 = xpart(1, i);
    x2 = xpart(2, i);
    x3 = xpart(3, i);
    x4 = xpart(4, i);
    x5 = xpart(5, i);
    h_r = 0.013;
    R_r = 0.0333;
    m = 0.063;
    %% state space
    f1 = x2;
    zm_in = x1 + x3;
    h_bar = (-zm_in + h_r)/(2.5*R_r);
    K_G = x4 + x5*(2/h_bar)^2;
    f2 = 9.81 - K_G/m*u_in(t);
    f3 = 0;
    f4 = 0;
    f5 = 0;
    f = [f1;f2;f3;f4;f5];
    xpart(:, i) = xpart(:, i) + f * Ts;
    xpartminus(:,i) = xpart(i);
    ypart = [1 0 0 0 0] * xpartminus(:,i);
    vhat = zm(t) - ypart;
    q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
end
% Normalize the likelihood of each a priori estimate.
    qsum = sum(q);
    for i = 1 : N
        q(i) = q(i) / qsum;
    end
%
 % Resample.
    for i = 1 : N
        u = rand; % uniform random number between 0 and 1
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + q(j);
            if qtempsum >= u
                xpart(i) = xpartminus(j);
                break;
            end
        end
    end
% The particle filter estimate is the mean of the particles.
xhatPart = mean(xpart,2);
%%
tnew=[tnew t];
xnew=[xnew x];
xhatnew=[xhatnew xhat];
P11new=[P11new P(1,1)];
%%
xhatpfnew=[xhatpfnew xhatPart];
end
%%
figure(1)
plot(tnew,xnew(1,:),'k',tnew,xhatpfnew(1,:),'b.');
%


% % figure(2)
% % plot(tnew,P11new,'k');

