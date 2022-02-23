%% Set parameters needed for the linear model

% Gravity constant
g = 9.81;       % m/s^2

% Wheel parameters
mw = 0.016;     % [kg] mass
rw = 0.021;     % [m] radius
Jw = mw*rw^2/2; % [kg*m^2] moment of inertia

% Body parameters
% Rough estimation (approximation by a cuboid), but well... it works :)
mb = 0.672;     % [kg] mass
Lc = 0.075;     % [m] distance to center of gravity
w = 0.11;       % [m] width
h = 0.09;       % [m] height
Jb = mb*w*h/8 + mb*Lc^2; % [kg*m^2] moment of inertia at wheel joint

% Motor parameters
Kb = 0.468;     % [V*s/rad] Back EMF constant
Kt = 0.317;     % [N*m/A] Torque constant
fm = 0.0022;    % Friction coefficient between body and motor
Jm = 10^-5;     % [kg*m^2] Motor moment of inertia
Rm = 6.69;      % [Ohm]
c = 9/100;      % [V] Coefficient for conversion PWM % to voltage

% Angle conversion
deg2rad = pi/180;

%% Set equations of motion
% M*x''+P*x'+Q*x = h    ': derivative
% x = [phi; alpha]

% Mass matrix
M = [rw^2*(mb+mw)+Jw+Jm,    mb*Lc*rw-Jm
     mb*Lc*rw-Jm,           mb*Lc^2+Jb+Jm]*deg2rad; 

% Matrix of velocity dependent forces
P = [Kt*Kb/Rm+fm,    -(Kt*Kb/Rm+fm)
     -(Kt*Kb/Rm+fm), Kt*Kb/Rm+fm]*deg2rad;

% Matrix of position dependent forces
Q = [0,  0
     0,  -mb*Lc*g]*deg2rad;

% Excitation forces
h = [Kt/Rm; -Kt/Rm]*c;

%% Derive state space representation
E = [eye(2),     zeros(2,2);
     zeros(2,2), M];

A = [zeros(2,2), eye(2)
     -Q,         -P];
 
B = [0; 0; h];
 
A = inv(E)*A;
B = inv(E)*B;

C = [1 0 0 0;
     0 1 0 0];
D = 0;

model = ss(A,B,C,D);

%% Design LQR controller
Q = [1 0 0 0;
     0 100 0 0;
     0 0 1 0;
     0 0 0 1];
 
% Since input to the LEGO Mindstorms is limited to [-100, 100], punish 
% input value more, such that they stay in the given limits.
R = 150;        

K = lqr(A,B,Q,R)

ctrl_sys = ss(A-B*K, B, C, D);
% step(ctrl_sys);