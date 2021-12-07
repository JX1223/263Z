% Clear
clear;
clc;
close all;
fs = 14;
%% Read geometry file to get x and y coordinates
%load("0012_199.mat"); af_name = "NACA 0012";
%load("4412_199.mat"); af_name = "NACA 4412";
load("4515_199.mat"); af_name = "NACA 4515";
X = coords(:,1);
Y = coords(:,2);

%% Loop through angle of attacks, -16 to 16
angles = -16:0.5:16;

x_sep_upper_all = zeros(1, numel(angles));
x_sep_lower_all = zeros(1, numel(angles));
CL_all = zeros(1, numel(angles));
CD_all = zeros(1, numel(angles));

% Create figure for plots
f1 = figure;
f1.Position(3:4) = [1500 500];
tl = tiledlayout(1,3, 'TileSpacing', 'tight');

nexttile;
hold on;
% Loop
for a = 1:numel(angles)
    %% Angle of attack
    alpha = deg2rad(angles(a));
    %% Apply vortex panel method
    [Cp, V, gamma, S, x, y, theta, m] = vortex_panel(X, Y, alpha);
    
    %% Calculate separation with Thwaites
    [x_sep_upper, x_sep_lower, LE_point] = get_separation(V, S, X, Y, x, y);
    
    % Store in array
    x_sep_upper_all(a) = x_sep_upper;
    x_sep_lower_all(a) = x_sep_lower;
    
    %% Calculate CP and plot
    if mod(angles(a), 8) == 0
        plot(X(floor(end/2):-1:1),-Cp(end/2:-1:1)+Cp(end/2+1:end),'LineWidth', 1, 'DisplayName', sprintf("\\alpha = %i\\circ", angles(a)));
    end
    
    %% Calculate CL
    % Get normal vectors
    nX = cos(theta + pi/2);
    nY = sin(theta + pi/2);
    CFx = sum(Cp.*(-nX).* S);
    CFy = sum(Cp.*(-nY).* S);
    R = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)]; % Rotation matrix
    CDCL = R * [CFx; CFy];
    CD_all(a) = CDCL(1); % Store coefficient of drag
    CL_all(a) = CDCL(2); % Store coefficient of lift
end
legend();
set(gca, 'YDir','reverse')
xlabel("x/c");
ylabel("\DeltaC_p")
grid on;
title("Pressure Distribution", 'FontSize', fs);
hold off;

%% Plot CL and CD
% Can get analytical CL from thin airfoil theory
name_char = char(af_name);
if name_char(6) == '0' && name_char(7) == '0'
    CL_analytical = 2*pi*deg2rad(angles);
    a_L0 = 0;
else
    D1 = str2double(name_char(6));
    D2 = str2double(name_char(7));
    p = D2/10;
    zm = D1/100;
    theta_m = acos(1 - 2*p);
    A0 = deg2rad(angles) - zm/(pi*(p^2)) * ((2*p - 1) * theta_m + sin(theta_m)) - (zm)/(pi*(1-p)^2)*((2*p-1)*(pi - theta_m) - sin(theta_m));
    A1 = (2*zm/(pi*p^2))*( (2*p - 1)*sin(theta_m) + (1/4)*(2*theta_m + sin(2*theta_m))) + (2 *zm)/(pi*(1-p)^2)*(pi/2 - (2*p - 1)*sin(theta_m) - 1/4 * (2*theta_m + sin(2*theta_m)));
    CL_analytical = 2*pi*(A0 + (1/2)*A1);
    a_L0 =  zm/(pi*(p^2)) * ((2*p - 1) * theta_m + sin(theta_m)) + (zm)/(pi*(1-p)^2)*((2*p-1)*(pi - theta_m) - sin(theta_m))-(1/2)*A1; % Zero lift angle
end

% Plot actual vs analytical
nexttile;
colororder(gca, {'blue','red'})
yyaxis left;
plot(angles, CL_all, '-b', 'LineWidth', 1);
hold on;
plot(angles, CL_analytical, '--k', 'LineWidth', 1);
%plot(angles, 2*pi*(deg2rad(angles) - a_L0), '.r', 'LineWidth', 1);
scatter(rad2deg(a_L0),0,'o','filled','MarkerFaceColor','r');
ylabel("C_L");
hold off;

yyaxis right;
plot(angles, CD_all, '-r', 'LineWidth', 1);
ylabel("C_D");
legend("C_L","Analytical",sprintf("\\alpha_{L0} = %0.3f\\circ",rad2deg(a_L0)),"C_D",'location','southeast');
ylim([-0.0125, 0.0125])
grid on;
xlabel("\alpha");
% Set tick marks to degree symbol
xt=get(gca,'xtick');
for k=1:numel(xt)
xt1{k}=sprintf('%d°',xt(k));
end
set(gca,'xticklabel',xt1);
title("C_L & C_D vs \alpha", 'FontSize', fs);
hold off;



%% Plot upper separation points against angle of attack
% Find stall angle with 1D interp
ind = find(x_sep_upper_all <= 0.2, 1, 'first');
stall_angle = interp1([x_sep_upper_all(ind); x_sep_upper_all(ind - 1)], [angles(ind); angles(ind - 1)], 0.2);


nexttile;
plot(angles, x_sep_upper_all, 'b', 'LineWidth', 1);
hold on;
scatter(stall_angle,0.2,'o', 'filled','MarkerFaceColor', 'r');
grid on;
title("Upper Separation", 'FontSize', fs);
ylabel("x_{s,u}/c");
xlabel("\alpha");
% Set tick marks to degree symbol
xt=get(gca,'xtick');
for k=1:numel(xt)
xt1{k}=sprintf('%d°',xt(k));
end
set(gca,'xticklabel',xt1);
legend("Separation Dist (Normalized)", sprintf("Approx. Stall Angle: %.2f\\circ", stall_angle));
hold off;

title(tl, af_name, 'FontSize', fs)
%exportgraphics(f1, sprintf("%s.eps", af_name));

%%  Fit C_L(alpha)
coeffs = polyfit(deg2rad(angles), CL_all,1); % Linear fit of C_L to alpha

