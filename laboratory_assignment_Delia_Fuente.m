%% Code for the Laboratory Assignment
% Author: Delia Fuente Pascual
% Date: 18/06/2021
% MSc Future Power Networks 2021/2022

clear 
clc
close all

%% Parameters initialisation
L_l1 = 0.0262; % pu
L_l2 = 0.0262; % pu
R_l1 = 0.0103; % pu
R_l2 = 0.0358; % pu
V_base = 381; % Vrms
S_base = 10; % kVA
f_base = 50; % Hz
I_base = 10^3*S_base/(sqrt(3)*V_base);

%% Design a power flow calculation algorithm to calculate the operating point.

Q4 = 0; % pu
Q5 = 0; % pu
wo = 1; % pu
Vo = 1; % pu
Dw = 100; % pu
Dv = 20; % pu
P4 = 1; % pu
P5 = 1; % pu

syms V1 V2 V3 w th1 th2 th3;

% Mismatches
f = @(V1,V2,V3,w,th1,th2,th3)... 
[P4-(Dw*(wo-w))+(V1^2)*(real(1/(R_l1+j*w*L_l1)))+V1*V2*((real(-1/(R_l1+j*w*L_l1)))*cos(th1-th2)+(imag(-1/(R_l1+j*w*L_l1)))*sin(th1-th2)),...
Q4-(Dv*(Vo-V1))-(V1^2)*(imag(1/(R_l1+j*w*L_l1)))+V1*V2*((real(-1/(R_l1+j*w*L_l1)))*sin(th1-th2)-(imag(-1/(R_l1+j*w*L_l1)))*cos(th1-th2)),...
-(Dw*(wo-w))+(V2^2)*(real(1/(R_l1+j*w*L_l1))+real(1/(R_l2+j*w*L_l2)))+V2*V1*((real(-1/(R_l1+j*w*L_l1)))*cos(th2-th1)+(imag(-1/(R_l1+j*w*L_l1)))*sin(th2-th1))+V2*V3*((real(-1/(R_l2+j*w*L_l2)))*cos(th2-th3)+(imag(-1/(R_l2+j*w*L_l2)))*sin(th2-th3)),...
-(Dv*(Vo-V2))-(V2^2)*(imag(1/(R_l1+j*w*L_l1))+imag(1/(R_l2+j*w*L_l2)))+V2*V1*((real(-1/(R_l1+j*w*L_l1)))*sin(th2-th1)-(imag(-1/(R_l1+j*w*L_l1)))*cos(th2-th1))+V2*V3*((real(-1/(R_l2+j*w*L_l2)))*sin(th2-th3)-(imag(-1/(R_l2+j*w*L_l2)))*cos(th2-th3)),...
P5-(Dw*(wo-w))+(V3^2)*(real(1/(R_l2+j*w*L_l2)))+V3*V2*((real(-1/(R_l2+j*w*L_l2)))*cos(th3-th2)+(imag(-1/(R_l2+j*w*L_l2)))*sin(th3-th2)),...
Q5-(Dv*(Vo-V3))-(V3^2)*(imag(1/(R_l2+j*w*L_l2)))+V3*V2*((real(-1/(R_l2+j*w*L_l2)))*sin(th3-th2)-(imag(-1/(R_l2+j*w*L_l2)))*cos(th3-th2))]';

fp = @(x) [f(x(1),x(2),x(3),x(4),x(5),x(6),x(7))]';

% Initial value of the iteration
X0=[1,1,1,1,0,0,0];

% Solution for the mismatches
X = fsolve(fp,X0);

%% Results
% Voltages
V = [X(1),X(2),X(3)];
disp('V1 (V) = ');
disp(V(1)*V_base);
disp('V2 (V) = ');
disp(V(2)*V_base);
disp('V3 (V) = ');
disp(V(3)*V_base);
V1 = V(1)*exp(j*X(5));
V2 = V(2)*exp(j*X(6));
V3 = V(3)*exp(j*X(7));

% Frequency
w = X(4);
disp('w (Hz) = ');
disp(w*f_base);

% Admittance matrix
Zl1 = R_l1+j*w*L_l1;
Zl2 = R_l2+j*w*L_l2;
Y = [1/Zl1, -1/Zl1, 0;-1/Zl1, 1/Zl1 + 1/Zl2, -1/Zl2 ; 0, -1/Zl2, 1/Zl2];

% Currents
I = V*Y;
disp('I12 (A) = ');
I12 = round((V1-V2)/Zl1,8)*I_base;
disp(I12);
disp('I32 (A) = ');
I32 = round((V3-V2)/Zl2,8)*I_base;
disp(I32)

% Power
disp('P1 = P2 = P3 (W) = ');
P = Dw*(wo-X(4))*S_base*1000;
disp(P);
disp('Q1 (var) = ');
Q1 = Dv*(Vo-X(1))*S_base*1000;
disp(Q1);
disp('Q2 (var) = ');
Q2 = Dv*(Vo-X(2))*S_base*1000;
disp(Q2);
disp('Q3 (var) = ');
Q3 = Dv*(Vo-X(3))*S_base*1000;
disp(Q3);
