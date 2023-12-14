%% parameters:
%% Sample rate
DeltaT = 0.01; %Sample time
%% Earth Parameters
f = 1/298.2572236;
e = 0.0818191908426;
R0 = 6378137;
ee = e*e;
Omega_e = 7.29211651452621e-5;
%% Sin and Cos Lateral
CCLA = cos(Pc(1,1));
SCLA = sin(Pc(1,1));
SC2LA = sin(2*Pc(1,1));
%% Calculating Ground Radious
RG = R0*(1-f*SCLA*SCLA); %Gussian Radius
RN = R0*(1-ee)/(1-ee*SCLA*SCLA)^1.5;
RE = R0/(1-ee*SCLA*SCLA)^0.5;
%% Sin and Cos Euler Angles
SCR = sin(Ac(1,1));
CCR = cos(Ac(1,1));
SCP = sin(Ac(2,1));
CCP = cos(Ac(2,1));
TCP = SCP/CCP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Body rate at NED frame from body rate at inertila frame
win(1,1) = Omega_e*CCLA + Vc(2,1)/(RE+Pc(3,1));
win(2,1) = -1* Vc(1,1)/(RN+Pc(3,1));
win(3,1) = -1*Omega_e*SCLA - Vc(2,1)*SCLA/((RE+Pc(3,1))*CCLA);
%%
% wnb = wib -Cbn'*win;
%% Calculating new attitude rate from body rate at NED
% ADot(1,1) = wnb(1,1)+(wnb(2,1)*SCR + wnb(3,1)*CCR)*TCP;
% ADot(2,1) = wnb(2,1)*CCR -wnb(3,1)*SCR;
% ADot(3,1) = (wnb(2,1)*SCR+ wnb(3,1)*CCR)/CCP;
%% Calculating new attitude angles
% Ac = Ac +ADot*DelatT;
%% Calculating Cos transfer Matrix for trasferring  from body frame to navigation frame
% Cbn=CalcCbn(Ac);
%% new method %%
% Omega_IB = [0 -wib(3) wib(2);
%             wib(3) 0 -wib(1) ;
%             -wib(2) wib(1) 0 ];
% Omega_IN = [0 -win(3) win(2);
%             win(3) 0 -win(1) ;
%             -win(2) win(1) 0 ];
% dot_Cbn = Cbn*Omega_IB - Omega_IN*Cbn;
% Cbn = Cbn + dot_Cbn*DeltaT;
% Ac(1, 1) = atan2(Cbn(3, 2), Cbn(3,3));
% Ac(2, 1) = asin(-Cbn(3, 1));
% Ac(3, 1) = atan2(Cbn(2, 1), Cbn(1, 1));
%% new mwthod2 
q0 = cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2);
q1 = sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2);
q2 = cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2);
q3 = cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2);

Cbn = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
       2*(q1*q2-q0*q3), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3+q0*q1);
       2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), q0^2-q1^2-q2^2+q3^2];
wnb = wib - Cbn'*win;

Omega_NB = [0 -wnb(3) wnb(2);
            wnb(3) 0 -wnb(1) ;
            -wnb(2) wnb(1) 0 ];
dot_q = 1/2*[0, -wnb';...
            wnb, -Omega_NB];
q = q + dot_q*DeltaT;
phi = atan2(2*(q(1)*q(2)+q(3)*q(4)), 1-2*(q(2)^2+q(3)^2));
theta = asin(2*(q(1)*q(3)-q(2)*q(4)));
psi = atan2(2*(q(1)*q(4)+q(2)*q(3)), 1-2*(q(3)^2+q(4)^2));

Ac = [phi; theta; psi]; 
%% Velocity Update
%% Transferring Force from body frame to NED frame
Fn=Cbn*Fb;
%% Calculating Gravity vector
gh0 = 9.780318*(1+5.3024e-3*SCLA*SCLA - 5.9e-6*SC2LA);
gh = gh0/(1+Pc(3,1)/RG)^2;
Gn(1,1)=0;
Gn(2,1)=0;
Gn(3,1)=gh;
%% Calculating Vdot
VDot(1,1)= Vc(1,1)*Vc(3,1)/(RN+Pc(3,1))-Vc(2,1)*Vc(2,1)*SCLA/((RE+Pc(3,1))*CCLA)...
           -2*Omega_e*Vc(2,1)*SCLA;
VDot(2,1)= 2*Omega_e*Vc(1,1)*SCLA + Vc(1,1)*Vc(2,1)*SCLA/((RE+Pc(3,1))*CCLA)...
           +2*Omega_e*Vc(3,1)*CCLA + Vc(2,1)*Vc(3,1)/(RE+Pc(3,1));
VDot(3,1)= -2*Omega_e*Vc(2,1)*CCLA - Vc(2,1)*Vc(2,1)/(RE+Pc(3,1))...
           -Vc(1,1)*Vc(1,1)/(RN+Pc(3,1));
VDot = Fn+VDot+Gn;
%% Calculating new linear velocity:euler discrtization
Vc = Vc + VDot*DeltaT;
%% Position Update:
PDot(1,1) = Vc(1,1)/(RN+Pc(3,1));
PDot(2,1) = Vc(2,1)/((RE+Pc(3,1))*CCLA);
PDot(3,1) = -1 * Vc(3,1);
%% Calculating new geodetic position
Pc = Pc+PDot*DeltaT;
%%


       
       






