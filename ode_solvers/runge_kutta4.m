function [t, y] = runge_kutta4(f, t0, y0, h, t_end)
% RUNGE_KUTTA4  Solve the ODE y' = f(t,y) with the classical 4th-order
%               Runge–Kutta method (RK4).
%
%   [t, y] = runge_kutta4(f, t0, y0, h, t_end)
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
%     [t, y] = runge_kutta4(@(t,y) -y, 0, 1, 0.1, 2);
%     plot(t, y, t, exp(-t))

    n = ceil((t_end - t0) / h);
    t = zeros(n+1, 1);
    y = zeros(n+1, 1);

    t(1) = t0;
    y(1) = y0;

    for k = 1:n
        tk = t(k);
        yk = y(k);

        k1 = f(tk,           yk);
        k2 = f(tk + h/2,     yk + h/2 * k1);
        k3 = f(tk + h/2,     yk + h/2 * k2);
        k4 = f(tk + h,       yk + h   * k3);

        t(k+1) = tk + h;
        y(k+1) = yk + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end
