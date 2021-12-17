function analyze_blade(a, c, c_f, const, fs)
    fp = ['.' filesep 'figures' filesep];
    %% Calculated values
    % Calculate radius
    dr = (const.r_t - const.r_r)/(numel(a)); % Width of each blade element
    r = linspace(const.r_r, const.r_t, numel(a)+1)'; 
    r = r(1:end-1) + dr/2; % Distance to 'midpoint' of each blade element
    
    % Calculate relative velocities and angles
    V_t = const.omega * r;
    V_r = sqrt(V_t.^2 + const.V_h^2);
    phi = atan(const.V_h./V_t);
    
    % Calculate dynamic pressure
    q = (1/2)*const.rho.*V_r.^2;

    %% Calculate results
    c_l = 2*pi*(a - const.a_L0); % sectional coefficient of lift
    c_n = c_l .* cos(phi); % section thrust coefficient
    C_N = c_n .* dr; % coefficient of thrust per blade element
    T = C_N .* q .* c; % Thrust per blade element
    thrust_blade = sum(T); % Thrust per blade is sum of thrust per blade element
    %
    c_t = c_l .* sin(phi); % section force coefficient, tangential
    C_T = c_t .* dr; % coefficient of force tangential
    D = C_T .* q .* c; % Drag per blade element
    Q = D .* r; % Torque per blade element
    Q_blade = sum(Q); % Torque per blade is sum of torque per blade element
    % Viscous effects
    c_t_visc = 2 .* c_f .* cos(phi);
    C_T_visc = c_t_visc .* dr;
    D_skin = C_T_visc .* q .* c; % Drag per surface of blade element
    Q_skin = D_skin .* r; % Torque caused by skin friction drag
    Q_visc = sum(Q_skin);
    Q_tot = Q_blade + Q_visc;

    %% Solidity
    blade_area = sum(c.*dr);
    solidity = (2*blade_area)/const.A_R;

    %% Power
    Q_rotor = 2 * Q_tot; % [N/m]
    P_rotor = Q_rotor * const.omega; % [W]
    P_tot = 4 * P_rotor;

    %% Print output
    fprintf("Thrust per blade: %0.5f [n]\n", thrust_blade);
    fprintf("Total torque per blade: %0.5f [n/m]\n", Q_tot);
    fprintf("\t Torque from induced drag: %0.5f [n/m]\n", Q_blade);
    fprintf("\t Torque from viscous drag: %0.5f [n/m]\n", Q_visc);
    fprintf("Solidity: %0.5f\n", solidity);
    fprintf("Total power: %0.5f [W] %0.3f [hp]\n", P_tot, P_tot/745.7);

    %% Generate figures
    f1 = figure;
    f1.Position = [100   300   800   600];
    tl = tiledlayout(3, 1);
    nexttile;
    plot(r, c_n, 'LineWidth', 1);
    ylabel("c_n", 'FontSize', fs);
    grid on;
    set(gca,'FontSize',fs)
    xlim([0, 10]);

    nexttile;
    plot(r, c_t, "LineWidth", 1);
    ylabel("c_t", 'FontSize', fs);
    grid on;
    set(gca,'FontSize',fs)
    xlim([0, 10]);

    nexttile;
    plot(r, c_t_visc, "LineWidth", 1);
    ylabel("c_{t,visc}", 'FontSize', fs);
    grid on;
    xlim([0, 10]);
    xlabel(tl,"Blade Length [m]", 'FontSize', fs);
    set(gca,'FontSize',fs)
    exportgraphics(f1, strcat(fp,'coeffs.eps'));

    %% Separate figure of loading per section
    f2 = figure;
    plot(r, T, 'LineWidth', 1);
    grid on;
    ylabel("Thrust [N]", 'FontSize',fs)
    xlabel("Blade Length [m]", 'FontSize', fs);
    set(gca,'FontSize',fs)
    exportgraphics(f2, strcat(fp,'loading.eps'));
end