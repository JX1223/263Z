function [Cp, V, gamma, S, xcp, ycp, theta, m] = vortex_panel(X, Y, alpha)
    %% Get number of panels
    m = length(X) - 1; 

    %% Get vortex panel strengths at panel points
    % Get truncated X and Y vectors (without the m+1 panel point coordinate)
    Xj = X(1:end-1);
    Yj = Y(1:end-1);

    % Calculate control points (midpoint of panels)
    x = zeros(m,1);
    y = zeros(m,1);
    for i = 1:m
        x(i) = (X(i) + X(i+1))/2;
        y(i) = (Y(i) + Y(i+1))/2;
    end

    % Build coefficient building blocks without loop
    theta = atan2(diff(Y), diff(X));
    theta_i = repmat(theta,1,m);
    theta_j = repmat(theta',m,1); % (m x m) matrix where each element (i,j) is theta(j)
    theta_i_j = theta_i - theta_i'; % (m x m) matrix where element (i,j) is theta(i) - theta(j)
    theta_i_2j = theta_i - 2*theta_i'; % (m x m) matrix where element (i,j) is theta(i) - 2*theta(j)
    xi_Xj = repmat(x, 1, m) - repmat(Xj', m, 1); % (m x m) matrix where element (i,j) is x(i) - X(j)
    yi_Yj = repmat(y, 1, m) - repmat(Yj', m, 1);
    S = sqrt(diff(X).^2 + diff(Y).^2); % Panel lengths
    Sj = repmat(S', m, 1); % element (i,j) is S(j)

    % Calculate coefficients Cn1_ij and Cn2_ij from geometry
    A = -1.*xi_Xj.*cos(theta_j) - (yi_Yj).*sin(theta_j);
    B = xi_Xj.^2 + yi_Yj.^2;
    C = sin(theta_i_j);
    D = cos(theta_i_j);
    E = xi_Xj.*sin(theta_j) - (yi_Yj).*cos(theta_j);
    F = log(1 + (Sj.^2 + 2.*A.*Sj)./B);
    G = atan2(E.*Sj,B+A.*Sj);
    P = (xi_Xj).*sin(theta_i_2j) + (yi_Yj).*cos(theta_i_2j);
    Q = (xi_Xj).*cos(theta_i_2j) - (yi_Yj).*sin(theta_i_2j);

    % Normal velocity constants
    Cn2 = D + (1/2).*Q.*F./Sj - (A.*C + D.*E).*G./Sj;
    Cn1 = (1/2).*D.*F + C.*G - Cn2;
    Cn1(eye(m) == 1) = -1;
    Cn2(eye(m) == 1) = 1;

    % Tangential velocity constants
    Ct2 = C + (1/2).*P.*F./Sj + (A.*D - C.*E).*G./Sj;
    Ct1 = (1/2).*C.*F - D.*G - Ct2;
    Ct1(eye(m) == 1) = pi/2;
    Ct2(eye(m) == 1) = pi/2;

    %% Calculate Gamma
    % Build matrix H such that H * gamma = RHS
    H = zeros(m+1);
    H(1:m, 1) = Cn1(:,1);
    for j = 2:m
        H(1:m,j) = Cn1(:,j) + Cn2(:,j-1);
    end
    H(1:m, m+1) = Cn2(:,m);
    % Kutta condition
    H(m+1,1) = 1;
    H(m+1,m+1) = 1;

    % Build RHS
    RHS = zeros(m+1, 1);
    RHS(1:m) = sin(theta - alpha);

    % Get circulation density
    gamma = H\RHS;

    %% Get tangential velocity
    % Build matrix K
    K = zeros(m,m+1);
    K(:, 1) = Ct1(:,1);
    for j = 2:m
        K(:,j) = Ct1(:,j) + Ct2(:,j-1);
    end
    K(:, m+1) = Ct2(:,m);

    % Calculate tangential normalized velocity
    V = cos(theta - alpha) + K*gamma;

    %% Calculate pressure coefficient
    Cp = 1 - V.^2; % Pressure coefficient at control points
    
    %% Export control points
    xcp = x;
    ycp = y;
end