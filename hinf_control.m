clear; clc;
m = 0.030;
Ixx = 1.43e-5;
Iyy = 1.43e-5;
Izz = 2.89e-5;
l = 0.046;
Kt = 2.3e-08;
Kd = 7.8e-11;
d = Kd/Kt;
g = 9.81;

Ar = [0         (d*g*m)/Kd      1               0                       0; ...
    -(d*g*m)/Kd     0           0               1                      0; ...
    0               0    -Kd/Ixx    (d*g*(Ixx-Izz)*m)/(Ixx*Kd)           0; ...
    0               0 -(d*g*(Ixx-Izz)*m)/(Ixx*Kd) -Kd/Ixx,      0; ...
    0               0           0              0                      -Kd/Izz];


Br = [0      0       0; ...
      0      0       0; ...
    l/(2*Ixx) 0     -l/(2*d*Ixx);...
    0       1/Ixx   0; ...
    0       0       1/Izz];


C = zeros(8,5);
%C = zeros(3,3);
%C(3,3) = 1;
C(1,1) = 1;
C(2,2) = 1;
C(3,5) = 1;
C(4,1) = 1;
C(5,2) = 1;
C(6,3) = 1;
C(7,4) = 1; 
C(8,5) = 1;

%C(1,1) = 1;
%C(2,2) = 1; 
%C(3,5) = 1;
%C(4,1) = 1;
%C(5,2) = 1; 
%C(6,5) = 1;

%C(4,3) = 1; 
%C(5,4) = 1;

P=ss(Ar,Br,C,zeros(8,3));

nmeas = 5;
ncon = 3;
[K,CL,asd] = hinfsyn(P,nmeas,ncon);
save('K_matrix.mat','K');
save('CL_matrix.mat','CL');