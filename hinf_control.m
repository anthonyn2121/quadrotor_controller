clear; clc;
%% Quadrotor parameters
%m = 0.030;
%Ixx = 1.43e-5;
%Iyy = 1.43e-5;
%Izz = 2.89e-5;
%l = 0.046;
%Kt = 2.3e-08;
%Kd = 7.8e-11;
d = 2.4e-3;
g = 9.81;

m = 1.79;            % kg
Ixx = 1.335e-2;      % kg m^-2
Iyy = 1.335e-2;      % kg m^-2
Izz = 2.465e-2;      % kg m^-2
l = 0.18;            % m
Kt = 8.82;           % N rad^-2 s^2
Kd = 1.09e-2;        % N m rad^-2 s^2


%% State dynamics
Ar = [0,         (d*g*m)/Kd,      1,               0,                       0; ...
    -(d*g*m)/Kd,     0,           0,               1 ,                     0; ...
    0,               0,    -Kd/Ixx,    (d*g*(Ixx-Izz)*m)/(Ixx*Kd),           0; ...
    0,               0, -(d*g*(Ixx-Izz)*m)/(Ixx*Kd), -Kd/Ixx,      0; ...
    0,               0,           0,              0,                      -Kd/Izz];


Br = [0,      0,       0; ...
      0,      0,       0; ...
    l/(2*Ixx), 0,     -l/(2*d*Ixx);...
    0,       1/Ixx,   0; ...
    0,       0,       1/Izz];


C = zeros(3,5);
%C = zeros(8,5);
%C(3,3) = 1;
C(1,1) = 1;
C(2,2) = 1;
C(3,5) = 1;

%C(4,3) = 1;
%C(5,4) = 1;
%C(4,1) = 1;
%C(5,2) = 1;
%C(6,3) = 1;
%C(7,4) = 1; 
%C(8,5) = 1;

P=ss(Ar,Br,C,zeros(3,3));
Cmat=C;

%% Hinf controller
nmeas = 5;
ncon = 3;

% Weighting functions
W1 = tf(1000) * mkfilter(0.2,2,'rc')  * mkfilter(2,1,'rc')^-1;
W3 = tf(1000) * mkfilter(0.2,2,'rc')^-1;
bode(W1);
figure(2)
bode(W3)
Paug = augw(P, [], [], W1);

%[K,CL,asd] = hinfsyn(Paug,nmeas,ncon);
[K,CL,asd] = loopsyn(P,W1);

figure(3);
sigma(P);
figure(4);
sigma(CL);
figure(5);
x = randn(27,1);
lsim(CL,zeros(3,1e2+1),linspace(0,1,1e2+1),x);

%% Robust feedback linearization

C = @(x) cos(x);
S = @(x) sin(x);
T = @(x) tan(x);

alpha = zeros(1,3);
alpha(1) = -2*l*Izz*Kd^2*(Ixx*(-.25*C(x(2))^(-2)*Izz*(2*d^2*g^2*m^2*(S(2*x(1))+x(1)+C(2*x(1))*(-4*S(x(1))+x(1))-S(2*x(2))*x(2)*T(x(1)))+ ...
    2*d*g*m*Kd*(-3*x(4)-4*C(2*x(2))*S(x(1))*x(5)+C(x(1))^(-1)*(2*x(4)-C(3*x(1))*x(4)+2*C(2*x(1))*C(2*x(2))*x(4)+S(3*x(1))*x(5))+(2*S(2*x(2))*x(3)+x(5))*T(x(1)))+ ...
        Kd^2*(-6*x(4)*x(5)-C(3*x(1))^(-1)*(2*S(2*x(2))*x(3)*x(4)+2*C(3*x(1))*x(4)*x(5)+S(3*x(1))*(x(4)-x(5))*(x(4)+x(5)))+ ...
            (x(5)^2-3*x(4)^2+2*C(2*x(2))*x(4)^2)*T(x(1)))-C(x(1))^(-1)*Kd^3*x(5)*T(x(2)))+Izz*(Kd^2*(C(x(1))^(-1)*Kd*x(5)*T(x(2))-d*g*m*(x(2)+C(x(1))^(-1)*...
                (-1+S(x(1))*x(1))*T(x(1))))) + Izz*(-Kd^2*x(4)*x(5)+d*g*m*Kd*(-x(4)+C(x(1))*x(4)-S(x(1))*x(5)+x(3)*T(x(1))*T(x(2))+S(x(1))*(x(5)+x(4)*T(x(1)))*T(x(2))^2))+...
                    d^2*g^2*m^2*(x(1)-x(2)*T(x(1))*T(x(2)) + S(x(1))*(-1+T(x(2))^2))));

alpha(2) = Izz*Kd^2*(Ixx*(-Kd^3*x(5)*T(x(1))+ Izz*(d^2*g^2*m^2*(C(x(1))^(-1)*x(2)+(-2+C(x(1)))*T(x(2)))+Kd^2*(2*x(5)+x(4)*T(x(1))+ ...
    C(x(1))^(-1)*(S(x(1))*x(4)+C(x(1))*x(5))^2*T(x(2)))+2*d*g*m*Kd*(x(3)-C(x(1))^(-1)*x(3)+(S(x(1))*x(4)-x(5)+C(x(1))*x(5)-x(4)*T(x(1)))*T(x(2))))+ ...
        Izz*(C(x(1))^(-1)*Kd^2*(d*g*m*(S(x(1))*x(1))+S(x(1))*Kd*x(5))+Izz*(-Kd^2*x(3)*x(5)+d^2*g^2*m^2*(-C(x(1))^(-1)*x(2)+T(x(2)))+ ...
            d*g*m*Kd*((-1+C(x(1))^(-1)*x(3)+(x(5)+x(4)*T(x(1)))*T(x(2))))))));

beta = [
    1   -(2*T(x(1))*T(x(2)))/l -(2*C(x(1))^(-1)*Ixx*T(x(2)))/(l*Izz);
    0   S(C(x(1)))              (Ixx*T(x(1)))/Izz;
    0   0                               1];



%% LQR controller
Q=diag([1 1 10 10 10]);
R=diag([0.1 1 5]);
[Gain, S, poles] = lqr(Ar, Br, Q, R);
CLsys = ss(Ar,-Br*Gain,eye(5,5),zeros(5,5));
G = tf(CLsys);

