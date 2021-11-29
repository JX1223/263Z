function [x_sep_upper, x_sep_lower, LE_point] = get_separation(V, S, X, Y, x, y)
    % Calculate the x dist of separation point, as projected onto the chord
    % Determine location of leading edge panel point:
    LE_ind1 = find(round(V,6) > 0, 1); % Find first panel with positive velocity
    LE_ind2 = find(round(V,6) < 0, 1, 'last');
    Ue_upper = V(round(V,6)>0); % Upper velocities, defined at control points starting from leading edge point
    Ue_lower = -1*flip(V(round(V,6)<0));
    m = numel(x); % The number of panels is equal to the number of control points
    upper_panels = LE_ind1:m;
    lower_panels = LE_ind2:-1:1;
    
    % Construct xs_upper, distance from stagnation point
    if LE_ind1 - LE_ind2 > 1 % For the case that the stagnation point is located directly on a control point
        S_prev = S(LE_ind1 - 1);
    else
        S_prev = 0;
    end
    s_upper = zeros(numel(upper_panels),1);
    s_cumulative = 0; % Keeps track of cumulative distance so far
    for i = 1:numel(upper_panels)
        s_cumulative = s_cumulative + S(upper_panels(i))/2 + S_prev/2;
        s_upper(i) = s_cumulative;
        S_prev = S(upper_panels(i));
    end
    
    % Construct xs_lower, distance from stagnation point
    if LE_ind1 - LE_ind2 > 1
        S_prev = S(LE_ind1 - 1);
    else
        S_prev = 0;
    end
    s_lower = zeros(numel(lower_panels),1);
    s_cumulartive = 0;
    for i = 1:numel(lower_panels)
        s_cumulartive = s_cumulartive + S(lower_panels(i))/2 + S_prev/2;
        s_lower(i) = s_cumulartive;
        S_prev = S(lower_panels(i));
    end
    
    % Get separation distance from LE stagnation point with Thwaites
    sep_upper = thwaites(s_upper, Ue_upper);
    sep_lower = thwaites(s_lower, Ue_lower);

    % Need to convert this distance to x coord projected on chord
    if y(upper_panels(1)) >= 0
        x_sep_upper = sep_upper + x(upper_panels(1));
    else
        x_sep_upper = sep_upper - x(upper_panels(1));
    end
    
    if y(lower_panels(1)) <= 0
        x_sep_lower = sep_lower + x(lower_panels(1));
    else
        x_sep_lower = sep_lower - x(lower_panels(1));
    end

    % Get leading edge point coordinates
    if LE_ind1 - LE_ind2 > 1
        LE_point = [x(LE_ind1 - 1); y(LE_ind1 - 1)];
    else
        LE_point = [X(LE_ind1); Y(LE_ind1)];
    end
end