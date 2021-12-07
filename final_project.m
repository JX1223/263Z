%
clear;
clc;
close all;

%% Inputs
% Constants
% Air constants
rho = 1.225; % kg/m^3
nu = 1.467*10^(-5);
% Problem constants
omega = 15 * (2*pi/60); % rpm -> rad/s
r_r = 1; % m, root radius
r_t = 10; % m, tip radius
a_L0 = -0.08; % rad, zero lift angle, from vortex panel
a_stall = deg2rad(8.473318733082900); % rad, stall angle, from vortex panel
mass = 115;
g = 9.81;
T_R = (mass * g)/(4); % N, thrust required per rotor
T_B = T_R / 2; % Thrust required per blade
A_R = pi*r_t^2; % m^2, area of disk
V_h = sqrt(T_R/(2*rho*A_R)); % Hover induced velocity

% Adjustable parameters
m = 101; % number of blade elements + 1

%% Calculated values
% Calculate radius
dr = (r_t - r_r)/(m-1); % Width of each blade element
r = linspace(r_r, r_t, m)'; 
r = r(1:end-1) + dr/2; % Distance to 'midpoint' of each blade element

% Calculate relative velocities and angles
V_t = omega * r;
V_r = sqrt(V_t.^2 + V_h^2);
phi = atan(V_h./V_t);

% Calculate constants used in dQ, dT
dQC = (1/2)*rho.*V_r.^2 .* r .* dr;
dTC = (1/2)*rho.*V_r.^2 .* dr .* cos(phi);

%% Solve optimization problem
% Define problem and variables
prob = optimproblem("ObjectiveSense","minimize");
a = optimvar('a', m-1, 'LowerBound', a_L0, 'UpperBound', a_stall - deg2rad(1)); % Optimization variable for angle of attack
coeffs = optimvar('coeffs', 5); % 5 for cubic, 4 for exponential

c_func = @(x) coeffs(1)*(x - coeffs(5)).^3 + coeffs(2)*(x - coeffs(5)).^2 + coeffs(3)*(x - coeffs(5)) + coeffs(4);
%c_func = @(x) coeffs(1)*exp(coeffs(2)*(x - coeffs(4))) + coeffs(3);

%c_ub = @(x) ((0.8 - 2.5)/9 * x) + (2.5 - (0.8 - 2.5)/9);
%c_lb = @(x) ((0.5 - 1.5)/9 * x) + (1.5 - (0.5 - 1.5)/9);
%c_ub = 2.5*ones(size(r)); c_ub(end) = 0.8;
%c_lb = 0.5*ones(size(r)); c_lb(1) = 1.5;
%c = optimvar('c', m-1, 'LowerBound', c_lb, 'UpperBound', c_ub); % Optimization variable for chord length

% Define objective function
expr = optimexpr;
for i = 1:m-1
    C_L = 2*pi*(a(i) - a_L0); % Sectional coefficient of lift
    c = c_func(r(i));
    Re_L = (V_r(i) * c)/nu; % Reynolds number w.r.t. chord
    C_f = 1.33/sqrt(Re_L); % Avg. skin friction coefficient, laminar
    dQ = ((C_L * sin(phi(i)) + 2*C_f)*c*dQC(i)); % Torque per element
    expr = expr + dQ;
end
prob.Objective = expr;

% Define Thrust constraint
constr1 = sum(2*pi*(a - a_L0) .* c_func(r) .* dTC) == T_B;
prob.Constraints.constr1 = constr1;

% Define chord constraint
constr2 = c_func(r) <= 2.5;
constr3 = c_func(1) >= 1.5;
constr4 = c_func(10) <= 0.8;
constr5 = c_func(r) >= 0.5;
prob.Constraints.constr2 = constr2;
prob.Constraints.constr3 = constr3;
prob.Constraints.constr4 = constr4;
prob.Constraints.constr5 = constr5;

% Define initial solution

m = (0.15 - 0.5)/9;
y = @(x) m*x + 0.539;

%x0.a = y(ones(size(a)));
x0.a = 0.05*ones(size(a));
x0.coeffs = ones(size(coeffs));
%x0.coeffs(2) = 0
%x0.c = ones(size(c));

% Solve problem
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end
options = optimoptions('fmincon','UseParallel',true);
[sol, fval] = solve(prob, x0, 'Options',options);

%% Plot result
figure;
plot(r, sol.a);
title("Angle of Attack");

coeffs = sol.coeffs;
c_func = @(x) coeffs(1)*(x - coeffs(5)).^3 + coeffs(2)*(x - coeffs(5)).^2 + coeffs(3)*(x - coeffs(5)) + coeffs(4);
%c_func = @(x) coeffs(1)*exp(coeffs(2)*(x - coeffs(4))) + coeffs(3);
c = c_func(r);
figure;
plot(r, c);
title("Chord");
