glvs;
ts = 0.005;
load flightimuavpr.mat;
% imuplot(imu); insplot(avpr);
avp = inspure(imu(1:1800/ts,:),avpr(1,:)','f');
avpcmpplot(avpr, avp);