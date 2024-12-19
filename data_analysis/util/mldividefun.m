%% Linear/non-linear fitting model
function [yfit, amplitude, c] = mldividefun(fun, y, p)
    % First calculate the exponential without scaling, then fit the linear
    % scaling parameter.
    % Solve Y = XA, where y: (y1 ... yn) = (1 ... 1, x1 ... xn)(a0, a1),
    % matrix form of y = a0 + a1*x. In this case, the matrix X is exp(-x/tau),
    % while Y = y.
    try
        X = [ones(length(y), 1) fun(p)];
    catch
        X = [ones(length(y), 1) fun(p)'];
    end
    A = X\y; % Solve linear equation: yfit * amplitude = y.
    amplitude = A(2);
    c = A(1);
    yfit = amplitude*fun(p) + c;
end