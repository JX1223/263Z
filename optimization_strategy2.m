%
clear;
clc;
close all;

%rng(13);
%% Constants
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

%% Adjustable parameters
m = 100; % number of blade elements

%% Calculated values
% Calculate radius
dr = (r_t - r_r)/(m); % Width of each blade element
r = linspace(r_r, r_t, m+1)'; 
r = r(1:end-1) + dr/2; % Distance to 'midpoint' of each blade element

% Calculate relative velocities and angles
V_t = omega * r;
V_r = sqrt(V_t.^2 + V_h^2);
phi = atan(V_h./V_t);

% Calculate q constant
q = (1/2)*rho.*V_r.^2;

% %% First, optimize for chord function
% % Define problem and variables
% prob = optimproblem("ObjectiveSense","minimize");
% a = optimvar('a', m, 'LowerBound', a_L0, 'UpperBound', a_stall - deg2rad(1)); % Optimization variable for angle of attack
% d = optimvar('d', m', 'LowerBound', 0, 'UpperBound', sqrt(2.5));
% % order_c = 4;
% % cfs = optimvar('cfs', order_c + 1);
% % c_func = @(x) build_polynomial(cfs, x);
% cfs = optimvar('cfs', 2);
% c_func = @(x) cfs(1) * exp(cfs(2)*x) + 0.5;
% 
% % Define objective function
% expr = optimexpr;
% for i = 1:m-1
%     C_L = 2*pi*(a(i) - a_L0); % Sectional coefficient of lift
%     D = 2*((1.33*nu^(1/2)*d(i))/(V_r(i)^(1/2)))*r(i)*q(i); % Drag from skin friction, Laminar
%     dQ = ((C_L * sin(phi(i))*c_func(r(i))*r(i)*q(i))+D)*dr; % Torque per element
%     expr = expr + dQ;
% end
% prob.Objective = expr;
% 
% % Define Thrust constraint
% constr1 = sum(2*pi*(a - a_L0) .* cos(phi) .* q .* c_func(r) .* dr) == T_B;
% prob.Constraints.constr1 = constr1;
% 
% % Define d and c constraint
% constr2 = d.^2 >= c_func(r);
% prob.Constraints.constr2 = constr2;
% 
% 
% constr3 = c_func(r) <= 2.5;
% constr4 = c_func(r) >= 0.5;
% constr5 = c_func(1) >= 1.5;
% constr6 = c_func(1) <= 2.5;
% constr7 = c_func(10) <= 0.8;
% constr8 = c_func(10) >= 0.5;
% 
% prob.Constraints.constr3 = constr3;
% prob.Constraints.constr4 = constr4;
% prob.Constraints.constr5 = constr5;
% prob.Constraints.constr6 = constr6;
% prob.Constraints.constr7 = constr7;
% prob.Constraints.constr8 = constr8;
% 
% % Define a constraint
% constr9 = a >= a_L0;
% constr10 = a <= a_stall - deg2rad(1);
% prob.Constraints.constr9 = constr9;
% prob.Constraints.constr10 = constr10;
% 
% % Define additional c constraints
% constr11 = cfs(2) >= -6;
% constr12 = cfs(2) <= -0.1;
% constr13 = cfs(1) >= 0.1;
% prob.Constraints.constr11 = constr11;
% prob.Constraints.constr12 = constr12;
% prob.Constraints.constr13 = constr13;
% 
% % Define initial solution 
% x0.a = rand(size(a));
% x0.cfs = rand(size(cfs));
% x0.d = rand(size(d));
% 
% % Solve problem
% options = optimoptions('fmincon','ConstraintTolerance',1e-10);
% [sol,fval] = solve(prob, x0, 'Solver','fmincon', 'Options',options);
% 
% %c_func = @(x) build_polynomial(sol.cfs, x);
% c_func = @(x) sol.cfs(1) * exp(sol.cfs(2)*x) + 0.5;
% 
% figure;
% plot(r, c_func(r), 'LineWidth', 1);
% title("c");

%% Using optimized c values, optimize a function for alpha
%c = c_func(r);
c_func = @(x) ones(size(r));
%c_func = @(x) exp(-(2).*(x - 1)) + 0.5;
c = c_func(r);

% Define problem and variables
prob2 = optimproblem("ObjectiveSense","minimize");
order_a = 1;
afs = optimvar('afs', order_a + 1);
a_func = @(x) build_polynomial(afs, x);

% Define objective function
expr = optimexpr;
for i = 1:m
    C_L = 2*pi*(a_func(r(i)) - a_L0); % Sectional coefficient of lift
    Re_L = V_r(i)*c(i)/nu;
    C_f = 1.33/sqrt(Re_L);
    w_bending = 0;
    bending = (C_L*cos(phi(i))*q(i)*c(i)*dr)*r(i); % We also want to minimize bending??
    dQ = ((C_L*sin(phi(i))*c(i)*r(i)*q(i)) + 2*C_f*cos(phi(i))*q(i)*c(i)*r(i))*dr; % Torque per element
    expr = expr + dQ + (w_bending * bending);
end
prob2.Objective = expr;

% Define Thrust constraint
constr1 = sum(2*pi*(a_func(r) - a_L0) .* cos(phi) .* q .* c .* dr) >= T_B;
prob2.Constraints.constr1 = constr1;

% Define a constraint
constr2 = a_func(r) >= a_L0;
constr3 = a_func(r) <= a_stall - deg2rad(1);
prob2.Constraints.constr2 = constr2;
prob2.Constraints.constr3 = constr3;

% Define initial solution 
clear('x0');
x0.afs = rand(size(afs));

% Solve problem
options = optimoptions('fmincon');
[sol2,fval2,EXITFLAG,OUTPUT,LAMBDA] = solve(prob2, x0,'Solver','fmincon');

% Plot
a_func = @(x) build_polynomial(sol2.afs, x);
figure;
plot(r, a_func(r), 'LineWidth', 1);
ylim([a_L0, a_stall]);
title("a");

% %% Verify results
a = a_func(r);
%c = c_func(r);
[t_blade, Q_tot, Q_blade, Q_visc, solidity] = verify(a, c, m);
fprintf("Cost: %.3f\nVerified Thrust: %.3f\nVerified Torque: %.3f\n", fval2, t_blade, Q_tot);

%blade3D(a, c, m);

function output = build_polynomial(coeffs, x)
    output = 0;
    order = numel(coeffs) - 1;
    for i = 0:order
        output = output + coeffs(i+1)*x.^i;
    end
end