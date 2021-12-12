classdef Constants
    properties (Constant = true)
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
        T_R = (Constants.mass * Constants.g)/(4); % N, thrust required per rotor
        T_B = Constants.T_R / 2; % Thrust required per blade
        A_R = pi*Constants.r_t^2; % m^2, area of disk
        V_h = sqrt(Constants.T_R/(2*Constants.rho*Constants.A_R)); % Hover induced velocity
    end
end