function A = A_func(in1,u)
%A_func
%    A = A_func(IN1,U)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    16-Jan-2024 16:42:05

x1 = in1(1,:);
x3 = in1(3,:);
x5 = in1(5,:);
t2 = x1.*1.201201201201201e+1;
t3 = x3.*1.201201201201201e+1;
t4 = t2+t3-5.2e+1./3.33e+2;
t5 = 1.0./t4.^3;
t6 = t5.*u.*x5.*1.525334858668192e+3;
A = reshape([0.0,t6,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,0.0,0.0,u.*(-1.0e+3./6.3e+1),0.0,0.0,0.0,0.0,1.0./t4.^2.*u.*(-6.349206349206349e+1),0.0,0.0,0.0],[5,5]);
end
