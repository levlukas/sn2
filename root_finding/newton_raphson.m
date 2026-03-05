function root = newton_raphson(f, df, x0, tol, max_iter)
% NEWTON_RAPHSON  Find a root of f using the Newton–Raphson method.
%
%   root = newton_raphson(f, df, x0, tol, max_iter)
%
%   Inputs:
%     f        - function handle
%     df       - handle for the derivative f'
%     x0       - initial guess
%     tol      - convergence tolerance (default 1e-6)
%     max_iter - maximum number of iterations (default 100)
%
%   Output:
%     root - approximate root of f
%
%   Example:
%     root = newton_raphson(@(x) x^2 - 2, @(x) 2*x, 1.0)
%     % returns approximately 1.4142

    if nargin < 4, tol = 1e-6; end
    if nargin < 5, max_iter = 100; end

    x = x0;
    for k = 1:max_iter
        fx = f(x);
        dfx = df(x);

        if abs(dfx) < eps
            error('newton_raphson: derivative is near zero at x = %g.', x);
        end

        x_new = x - fx / dfx;

        if abs(x_new - x) < tol
            root = x_new;
            return;
        end

        x = x_new;
    end

    root = x;
    warning('newton_raphson: maximum iterations reached; result may not have converged.');
end
