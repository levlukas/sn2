function [t, y] = euler(f, t0, y0, h, t_end)
% EULER  Solve the ODE y' = f(t,y) with explicit Euler's method.
%
%   [t, y] = euler(f, t0, y0, h, t_end)
%
%   Inputs:
%     f     - function handle f(t, y) for the right-hand side
%     t0    - initial time
%     y0    - initial value y(t0)
%     h     - step size
%     t_end - end time
%
%   Outputs:
%     t - column vector of time points
%     y - column vector of solution values
%
%   Example:
%     % Solve y' = -y,  y(0) = 1  (exact: y = e^{-t})
%     [t, y] = euler(@(t,y) -y, 0, 1, 0.1, 2);
%     plot(t, y, t, exp(-t))

    n = ceil((t_end - t0) / h);
    t = zeros(n+1, 1);
    y = zeros(n+1, 1);

    t(1) = t0;
    y(1) = y0;

    for k = 1:n
        t(k+1) = t(k) + h;
        y(k+1) = y(k) + h * f(t(k), y(k));
    end
end
