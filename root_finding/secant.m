function root = secant(f, x0, x1, tol, max_iter)
% SECANT  Find a root of f using the secant method.
%
%   root = secant(f, x0, x1, tol, max_iter)
%
%   Inputs:
%     f        - function handle
%     x0, x1  - two initial guesses (need not bracket the root)
%     tol      - convergence tolerance (default 1e-6)
%     max_iter - maximum number of iterations (default 100)
%
%   Output:
%     root - approximate root of f
%
%   Example:
%     root = secant(@(x) cos(x) - x, 0, 1)
%     % returns approximately 0.7391

    if nargin < 4, tol = 1e-6; end
    if nargin < 5, max_iter = 100; end

    f0 = f(x0);
    f1 = f(x1);

    for k = 1:max_iter
        if abs(f1 - f0) < eps
            error('secant: denominator too small; method failed.');
        end

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if abs(x2 - x1) < tol
            root = x2;
            return;
        end

        x0 = x1;  f0 = f1;
        x1 = x2;  f1 = f(x2);
    end

    root = x1;
    warning('secant: maximum iterations reached; result may not have converged.');
end
