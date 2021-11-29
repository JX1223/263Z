function dydx = differentiate(x, y)
    % Differentiates values y wrt to x using the 3-point rule
    % Input: x, y of dimension n
    % Output: dydx of dimension n
    dydx = zeros(size(y));
    for i = 1:numel(y)
        if i == 1 % Use the left endpoint formula for the first point
            h1 = x(i+1) - x(i);
            h2 = x(i+2) - x(i+1);
            f0 = y(i);
            f1 = y(i+1);
            f2 = y(i+2);
            dydx(i) = -(2*h1+h2)/(h1*(h1+h2))*f0 + ((h1+h2)/(h1*h2))*f1 - (h1/((h1+h2)*h2))*f2;
        elseif i == numel(y) % Use the right endpoint formula for last point
            h2 = x(i) - x(i - 1);
            h1 = x(i - 1) - x(i - 2);
            f0 = y(i - 2);
            f1 = y(i - 1);
            f2 = y(i);
            dydx(i) = (h2/(h1*(h1+h2)))*f0 - ((h1+h2)/(h1*h2))*f1 + ((h1+2*h2)/(h2*(h1+h2)))*f2;
        else % Use the midpoint formula for interior points
            h1 = x(i) - x(i - 1);
            h2 = x(i+1) - x(i);
            f0 = y(i - 1);
            f1 = y(i);
            f2 = y(i + 1);
            dydx(i) = -(h2/(h1*(h1+h2)))*f0 - ((h1-h2)/(h1*h2))*f1 + (h1/(h2*(h1+h2)))*f2;
        end
    end
end