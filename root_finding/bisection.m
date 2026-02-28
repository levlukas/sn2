function root = bisection(f, a, b, tol, max_iter)
% BISECTION  Find a root of f in [a,b] using the bisection method.
%
%   root = bisection(f, a, b, tol, max_iter)
%
%   Inputs:
%     f        - function handle, e.g. @(x) x^2 - 2
%     a, b     - bracket endpoints with f(a)*f(b) < 0
%     tol      - convergence tolerance (default 1e-6)
%     max_iter - maximum number of iterations (default 100)
%
%   Output:
%     root - approximate root of f
%
%   Example:
%     root = bisection(@(x) x^3 - x - 2, 1, 2, 1e-8, 200)
%     % returns approximately 1.5214

    if nargin < 4, tol = 1e-6; end
    if nargin < 5, max_iter = 100; end

    if f(a) * f(b) > 0
        error('bisection: f(a) and f(b) must have opposite signs.');
    end

    for k = 1:max_iter
        c = (a + b) / 2;
        fc = f(c);

        if abs(fc) < tol || (b - a) / 2 < tol
            root = c;
            return;
        end

        if f(a) * fc < 0
            b = c;
        else
            a = c;
        end
    end

    root = (a + b) / 2;
    warning('bisection: maximum iterations reached; result may not have converged.');
end
