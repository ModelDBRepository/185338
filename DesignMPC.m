% DesignMPC.m
% This code requires the Multi-Parametric Toolbox
% (http://control.ee.ethz.ch/~mpt/).

clear all;

%% arm_parameters -----------------------------------------------
% Segment length
l1_ = 0.15;
l2_ = 0.21;

% Segment center of mass
s1_ = l1_/2;
s2_ = 0.11;

% Moment arms
J = [2.6 -1.3 0 0 0.7 -2.5;
   0 0 1.2 -1.7 1.6 -1.1]/100;

% Hand position
x0 = 0; y0 = 0.24;

%% Inverse kinematics--------------------------------------------
l0_ = sqrt(x0^2 + y0^2);

cos_theta2 = -((x0^2+y0^2)-(l1_^2+l2_^2))/(2*l1_*l2_);
sin_theta2 = sqrt((2*l1_*l2_)^2-(l0_^2 - (l1_^2+l2_^2))^2)/(2*l1_*l2_);

theta2 = atan2(real(sin_theta2),real((-cos_theta2)));

kc = l1_ + l2_*cos(theta2);
ks = l2_*sin(theta2);

cos_theta1 = (kc*x0+ks*y0)/(kc^2+ks^2);
sin_theta1 = (-ks*x0+kc*y0)/(kc^2+ks^2);

theta1 = atan2(sin_theta1,cos_theta1);
theta = [theta1 theta2];

% Jacobian matrix
Jacob =[-l1_*sin(theta1)-l2_*sin(theta1+theta2) -l2_*sin(theta1+theta2);
    l1_*cos(theta1)+l2_*cos(theta1+theta2) l2_*cos(theta1+theta2)];

%% Dynamics model -----------------------------------------------
dt = 0.01;
delta = 0.04;

A = zeros(8,8);
A(1,1) = 1 - dt/delta;
A(2,2) = 1 - dt/delta;
A(1:2,3:8) = dt/delta*J;
A(3:8,3:8) = (1 - dt/delta)*eye(6);

B = zeros(8,6);
B(3:8,1:6) =  dt/delta*eye(6);

C = zeros(4,8);
C(1:2,1:2) = eye(2);
C(3:4,1:2) = inv(Jacob');
D = zeros(4,6);

szX = size(A,1);
szU = size(B,2);
szY = size(C,1);

% Cost weight
cost.wt = 1e+3; % Torque
cost.r = 1; % Efort

q = zeros(2,szY);
q(1:2,1:2) = cost.wt*eye(2);
Q = q'*q;
R = cost.r * eye(szU);

model = LTISystem('A', A, 'B', B, 'C', C, 'D', D);
model.u.min = zeros(6,1);
model.u.max = [Inf; Inf; Inf; Inf; Inf; Inf];
model.y.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
model.y.with('reference');
model.y.reference = 'free';

%% Model predictive controller ----------------------------------
N = 3;
mpc = MPCController(model, N);
ctrl = mpc.toExplicit();

%% Optimizer ----------------------------------------------------
rNum = ctrl.optimizer.Num;
for i=1:rNum
    F{i} = ctrl.optimizer.Set(i).Functions('primal').F;
    g{i} = ctrl.optimizer.Set(i).Functions('primal').g;
    AA{i} = ctrl.optimizer.Set(i).A;
    bb{i} = ctrl.optimizer.Set(i).b;
end

save data;

