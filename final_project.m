%% Initialize
clear;
clc;
close all;
rng(12);
%% Get Constants
const = Constants;

%% Adjustable parameters
m = 100; % number of blade elements
fs = 14;

%% Calculated values
% Calculate radius
dr = (const.r_t - const.r_r)/(m); % Width of each blade element
r = linspace(const.r_r, const.r_t, m+1)'; 
r = r(1:end-1) + dr/2; % Distance to 'midpoint' of each blade element

% Calculate relative velocities and angles
V_t = const.omega * r;
V_r = sqrt(V_t.^2 + const.V_h^2);
phi = atan(const.V_h./V_t);

% Calculate dynamic pressure q
q = (1/2)*const.rho.*V_r.^2;

%% Optimize for angle of attack, alpha
% Define an arbitrary chord function
%c_func = @(x) ((0.8 - 2.5)/9).*x + (2.5 - (0.8 - 2.5)/9);
%c_func = @(x) ((0.5 - 1.5)/9).*x + (1.5 - (0.5 - 1.5)/9);
c_func = @(x) exp(-(2).*(x - 1)) + 0.5;
c = c_func(r);

% Get local drag coefficients
c_f = drag_coeffs(V_r, c, const);

% Define problem and variables
prob = optimproblem("ObjectiveSense","minimize");
order_a = 2; % Should be <= 3
afs = optimvar('afs', order_a + 1);
w_bending = 10;
a_func = @(x) build_polynomial(afs, x);

% Define objective function
expr = optimexpr;
for i = 1:m
    c_l = 2*pi*(a_func(r(i)) - const.a_L0); % Sectional coefficient of lift
    c_n = c_l .* cos(phi(i)); % Sectional coefficient of force, normal
    c_t = c_l .* sin(phi(i)); % Sectional coefficient of force, tangential
    C_T = c_t .* dr; % Coefficient of force, tangential
    C_N = c_n .* dr; % Coefficient of force, normal
    D = C_T .* q(i) .* c(i); % Induced drag
    c_t_visc = 2 * c_f(i) .* cos(phi(i)); % Sectional coeff of visc force, tan
    C_T_visc = c_t_visc .* dr; % Coefficient of viscous force, tangential
    D_visc = (C_T_visc .* q(i) .* c(i)); % Drag due to viscous force
    dQ = (D + D_visc).* r(i);
    T = C_N .* q(i) .* c(i);
    dM_bending = T .* r(i);
    
    expr = expr + dQ + (w_bending * dM_bending);
end
prob.Objective = expr;

% Define Thrust constraint
constr1 = sum(2*pi*(a_func(r) - const.a_L0) .* cos(phi) .* q .* c .* dr) >= const.T_B;
prob.Constraints.constr1 = constr1;

% Define angle of attack constraints
constr2 = a_func(r) >= const.a_L0;
constr3 = a_func(r) <= const.a_stall - deg2rad(1);
constr4 = a_func(1) >= const.a_L0;
constr5 = a_func(1) <= const.a_stall - deg2rad(1);
constr6 = a_func(10) >= const.a_L0;
constr7 = a_func(10) <= const.a_stall - deg2rad(1);
prob.Constraints.constr2 = constr2;
prob.Constraints.constr3 = constr3;
prob.Constraints.constr4 = constr4;
prob.Constraints.constr5 = constr5;
prob.Constraints.constr6 = constr6;
prob.Constraints.constr7 = constr7;

% Define initial solution 
clear('x0');
x0.afs = rand(size(afs));

% Solve problem
options = optimoptions('fmincon','Display','final-detailed', 'ConstraintTolerance',1e-10);
[sol,fval,EXITFLAG,OUTPUT,LAMBDA] = solve(prob, x0,'Solver','fmincon','Options',options);
a_func = @(x) build_polynomial(sol.afs, x);

%% Analyze results
a = a_func(r);
analyze_blade(a, c, c_f, const, fs);

%% Plot
r_plot = linspace(1, 10, 101);
V_t_plot = const.omega * r_plot;
phi_plot = atan(const.V_h./V_t_plot);

% Plot angle of attack
f1 = figure;
f1.Position = [100   300   1000   400];
tl1 = tiledlayout(1,2,'Padding','compact');
nexttile;
plot(r_plot, rad2deg(a_func(r_plot)), 'LineWidth', 1);
hold on;
plot(r_plot, rad2deg(const.a_stall)*ones(size(r_plot)), '--r', 'LineWidth', 1)
plot(r_plot, rad2deg(const.a_L0)*ones(size(r_plot)), '-.b', 'LineWidth', 1)
legend("\alpha","\alpha_{stall}","\alpha_{L0}",'Location','east','FontSize', fs);
hold off;
ylabel("\alpha",'FontSize', fs)
xlabel("Blade Radius [m]", 'FontSize', fs)
ylim([-5, 9]);
xlim([0, 10]);

yticks(-5:1:9);
yt=get(gca,'ytick');
yt1 = cell(numel(yt),1);
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);

grid on;
title("Angle of Attack \alpha", 'FontSize', fs);

% Plot section angle
nexttile;
plot(r_plot, rad2deg(a_func(r_plot) + phi_plot), 'LineWidth', 1);
hold on;
plot(r_plot, rad2deg(phi_plot), '--r', 'LineWidth', 1);
hold off;
legend("\beta","\phi", 'location','east', 'FontSize', fs);
ylabel("\beta, \phi",'FontSize', fs)
xlabel("Blade Radius [m]", 'FontSize', fs)
ylim([-5, 25]);
xlim([0, 10]);

yt=get(gca,'ytick');
yt1 = cell(numel(yt),1);
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);

grid on;
title("Section Angle \beta and Induced Angle \phi", 'FontSize', fs);

% Plot chord variation
f2 = figure;
f2.Position = [100   300   500   400];
plot(r_plot, c_func(r_plot), 'LineWidth', 1);
grid on;
ylabel("c [m]");
xlim([0, 10]);
xlabel("Blade Radius [m]", 'FontSize', fs)
title("Chord",'FontSize', fs);

blade3D(a, c, const);

function output = build_polynomial(coeffs, x)
    output = 0;
    order = numel(coeffs) - 1;
    for i = 0:order
        output = output + coeffs(i+1)*x.^i;
    end
end

function c_f = drag_coeffs(V_r, c, const)
    Re_L = (V_r .* c)./const.nu;
    c_f = zeros(size(c));
    c_f(Re_L <= 5 * 10^5) = 1.328 .* Re_L(Re_L <= 5*10^5).^(-1/2);
    c_f(Re_L > 5 * 10^5) = 0.074 .* Re_L(Re_L > 5 * 10^5).^(-1/5); 
end