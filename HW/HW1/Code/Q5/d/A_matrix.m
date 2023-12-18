function A = A_matrix(in1,in2,g,R_m,R_p,R_g,Omega_e)
%A_matrix
%    A = A_matrix(IN1,IN2,G,R_m,R_p,R_g,Omega_e)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    18-Dec-2023 22:52:58

L = in1(7,:);
a_x = in2(4,:);
a_y = in2(5,:);
a_z = in2(6,:);
h = in1(9,:);
phi = in1(1,:);
psi = in1(3,:);
theta = in1(2,:);
v_d = in1(6,:);
v_e = in1(5,:);
v_n = in1(4,:);
w_y = in2(2,:);
w_z = in2(3,:);
t2 = cos(L);
t3 = sin(L);
t4 = tan(L);
t5 = cos(phi);
t6 = cos(psi);
t7 = cos(theta);
t8 = sin(phi);
t9 = sin(psi);
t10 = sin(theta);
t11 = tan(theta);
t12 = R_m+h;
t13 = R_p+h;
t14 = v_e.^2;
t17 = 1.0./R_g;
t15 = t4.^2;
t16 = t11.^2;
t18 = t5.*t6;
t19 = 1.0./t2;
t20 = t5.*t9;
t21 = t6.*t8;
t22 = t8.*t9;
t23 = 1.0./t7;
t25 = Omega_e.*t2.*2.0;
t26 = Omega_e.*t3.*2.0;
t27 = 1.0./t12;
t29 = 1.0./t13;
t24 = t23.^2;
t28 = t27.^2;
t30 = t29.^2;
t31 = t15+1.0;
t32 = t16+1.0;
t33 = t10.*t22;
t34 = t10.*t18;
t35 = t10.*t20;
t36 = t10.*t21;
t37 = -t35;
t38 = -t36;
t39 = t18+t33;
t40 = t22+t34;
t41 = t27.*t31.*v_e;
t42 = t20+t38;
t43 = t21+t37;
t44 = t25+t41;
mt1 = [t5.*t11.*w_y-t8.*t11.*w_z,-t8.*w_y-t5.*w_z,t5.*t23.*w_y-t8.*t23.*w_z,a_y.*t40+a_z.*t42,-a_y.*t43-a_z.*t39,a_y.*t5.*t7-a_z.*t7.*t8,0.0,0.0,0.0,t8.*t32.*w_y+t5.*t32.*w_z,0.0,t8.*t10.*t24.*w_y+t5.*t10.*t24.*w_z,-a_x.*t6.*t10+a_y.*t7.*t21+a_z.*t7.*t18,-a_x.*t9.*t10+a_y.*t7.*t22+a_z.*t7.*t20,-a_x.*t7-a_y.*t8.*t10-a_z.*t5.*t10,0.0,0.0,0.0,0.0,0.0,0.0,-a_y.*t39+a_z.*t43-a_x.*t7.*t9,-a_y.*t42+a_z.*t40+a_x.*t6.*t7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t27.*v_d,t26+t4.*t27.*v_e,t27.*v_n.*-2.0,t27,0.0,0.0,0.0,0.0,0.0];
mt2 = [-t26-t4.*t27.*v_e.*2.0,t29.*v_d+t4.*t27.*v_n,-t25-t29.*v_e.*2.0,0.0,t19.*t29,0.0,0.0,0.0,0.0,t27.*v_n,t25+t29.*v_e,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,-t44.*v_e,t44.*v_n-Omega_e.*t3.*v_d.*2.0,t26.*v_e,0.0,t3.*t19.^2.*t29.*v_e,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t14.*t28-t28.*v_d.*v_n,-t30.*v_d.*v_e-t4.*t28.*v_e.*v_n,t14.*t30+t28.*v_n.^2+g.*t17.*(h.*t17+1.0).*2.0,-t28.*v_n,-t19.*t30.*v_e,0.0];
A = reshape([mt1,mt2],9,9);
end
