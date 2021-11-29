function x_bar = thwaites(x,Ue)
%% Constants
threshold = -0.09;

%% Calculate Integral
int_Ue = cumtrapz(x, Ue.^5);

%% Calculate derivative
dUedx = differentiate(x, Ue); % Differentiate with the 3 point method

%% Calculate K
K = (0.45./(Ue.^6)) .* dUedx .* int_Ue;
kx = x;
% Get x_bar, where k first crosses or equals threshold
if K(1) == threshold
    x_bar = K(1);
else
    if K(1) < threshold
        ind2 = find(K >= threshold, 1);
        if isempty(ind2) % K does not cross or equal threshold
            x_bar = [];
            return;
        end
        ind1 = ind2 - 1;
    elseif K(1) > threshold
        ind2 = find(K <= threshold, 1);
        if isempty(ind2) % K does not cross or equal threshold
            x_bar = [];
            return;
        end
        ind1 = ind2 - 1;
    end
    slope = (K(ind2) - K(ind1))/(kx(ind2)-kx(ind1));
    x_bar = (threshold - K(ind1))/slope + kx(ind1);
end
end