clear; clc;
m = 0.030;
Ixx = 1.43e-5;
Iyy = 1.43e-5;
Izz = 2.89e-5;
l = 0.046;
Kt = 2.3e-08;
Kd = 7.8e-11;
d = 2.4e-3;
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


C = zeros(3,5);
%C = zeros(3,3);
%C(3,3) = 1;
C(1,1) = 1;
C(2,2) = 1;
C(3,5) = 1;
%C(4,1) = 1;
%C(5,2) = 1;
%C(6,3) = 1;
%C(7,4) = 1; 
%C(8,5) = 1;

P=ss(Ar,Br,C,zeros(3,3));

nmeas = 3;
ncon = 3;

% Weighting functions
W1 = tf(1000) * mkfilter(0.2,2,'rc')  * mkfilter(2,1,'rc')^-1;
%W3 = tf(1000) * mkfilter(0.2,2,'rc')^-1;
Paug = augw(P, W1, [], []);

[K,CL,asd] = hinfsyn(Paug,nmeas,ncon);
bode(W1);
figure(2);
sigma(P);
figure(3);
lsim(CL,zeros(3,10001),linspace(0,20,10001),randn(14,1));

u = [-424941; 180967; -63065];

Q=diag([1 1 10 10 10]);
R=diag([0.1 1 5]);
[Gain, S, poles] = lqr(Ar, Br, Q, R);
CLsys = ss(Ar,-Br*Gain,eye(5,5),zeros(5,5));
G = tf(CLsys);
Gain

