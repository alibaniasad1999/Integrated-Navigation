function F = F_func(in1,in2)
%F_func
%    F = F_func(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    18-Jan-2024 10:10:31

b = in1(4,:);
p = in2(1,:);
phi = in1(1,:);
q = in2(2,:);
r = in2(3,:);
theta = in1(2,:);
t2 = cos(phi);
t3 = sin(phi);
t4 = b+q;
t5 = b+r;
t6 = t2.*t5;
t7 = t3.*t4;
t8 = t6+t7;
F = [b+p+t8.*tan(theta).^2;t2.*t4-t3.*t5;t8./cos(theta);b.*(-1.0./3.82e+2)];
end