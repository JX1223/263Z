function [thrust_blade, Q_tot, Q_blade, Q_visc, solidity] = verify(a, c, m)
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
    
    %% Calculated values
    % Calculate radius
    dr = (r_t - r_r)/(numel(a)); % Width of each blade element
    r = linspace(r_r, r_t, numel(a)+1)'; 
    r = r(1:end-1) + dr/2; % Distance to 'midpoint' of each blade element
    
    % Calculate relative velocities and angles
    V_t = omega * r;
    V_r = sqrt(V_t.^2 + V_h^2);
    phi = atan(V_h./V_t);
    
    % Calculate dynamic pressure
    q = (1/2)*rho.*V_r.^2;

    %% Calculate results
    c_l = 2*pi*(a - a_L0); % sectional coefficient of lift
    C_L = c_l .* dr; % coefficient of lift per blade element
    C_T = C_L .* cos(phi); % coefficient of thrust per blade element
    T = C_T .* q .* c; % Thrust per blade element
    thrust_blade = sum(T); % Thrust per blade is sum of thrust per blade element
    %
    C_D = C_L .* sin(phi); % coefficient of drag per blade element
    D = C_D .* q .* c; % Drag per blade element
    Q = D .* r; % Torque per blade element
    Q_blade = sum(Q); % Torque per blade is sum of torque per blade element
    % Viscous effects
    Re_L = (V_r.*c)./nu;
    C_f = 1.33./sqrt(Re_L); % Coefficient of skin friction
    D_skin = C_f .* cos(phi) .* q .* c .* dr; % Drag per surface of blade element
    Q_skin = 2 * D_skin .* r; % Torque caused by skin friction drag
    Q_visc = sum(Q_skin);
    Q_tot = Q_blade + Q_visc;

    %% Solidity
    blade_area = sum(c.*dr);
    solidity = (2*blade_area)/A_R;
end