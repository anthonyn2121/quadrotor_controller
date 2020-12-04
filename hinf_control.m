clear; clc;
%% Quadrotor parameters

d = 2.4e-3;
g = 9.81;
m = 1.79;            % kg
Ixx = 1.335e-2;      % kg m^-2
Iyy = 1.335e-2;      % kg m^-2
Izz = 2.465e-2;      % kg m^-2
l = 0.18;            % m
Kt = 8.82;           % N rad^-2 s^2
Kd = 1.09e-2;        % N m rad^-2 s^2
f_max = 1e-1;        % N
f_min = 0;
fdot_max = f_max * 10;     % N s^-1

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

P=ss(Ar,Br,C,zeros(3,3));
Cmat=C;

%% Hinf controller

% Weighting functions
W1 = tf(1000) * mkfilter(0.2,2,'rc')  * mkfilter(2,1,'rc')^-1;
bode(W1);

[K,CL,gam,INFO] = loopsyn(P,W1);

figure(2);
sigma(P);
figure(3);
sigma(CL);
figure(4);
x = randn(27,1);
lsim(CL,zeros(3,1e2+1),linspace(0,5,1e2+1),x);

%% Reduced Order controller
size(K)
figure(5);
hsv = hankelsv(K);
semilogy(hsv,'*--')
grid
title('Hankel singular values of K')
xlabel('Order')

figure(6);
Kr = reduce(K,8);
size(Kr)
sigma(K,'b',K-Kr,'r-.',{1e-4,1e6})
legend('K','error K-Kr')

figure(7);
Tr = feedback(P*Kr,eye(3));
T = feedback(P*K,eye(3));
step(T,'b',Tr,'r-.',8)
title('Responses to step commands for alpha and theta');
legend('K','Kr')
%% Robust feedback linearization
%v = K.D * x_hat + K.C * 
u = linearize(x,v)

%% Thrust from inputs u

F = 1/4 * [1, -2/l, 1/d; 1, 2/l, 1/d; 2, 0, -2/d] * u

set(gca,'FontSize', 20)
legend('Target loop shape','Singular Values of the Plant','Singular values of the shaped plant')
xlabel('Frequency ','FontSize',20)
ylabel('Singular Values','FontSize',20)
title('Singular Values and Target loop Shape')
title('Singular Values and Target loop Shape','FontSize',24)

