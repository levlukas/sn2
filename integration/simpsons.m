function I = simpsons(f, a, b, n)
% SIMPSONS  Approximate the integral of f over [a,b] using the
%           composite Simpson's 1/3 rule.
%
%   I = simpsons(f, a, b, n)
%
%   Inputs:
%     f    - function handle
%     a, b - integration limits
%     n    - number of sub-intervals (must be even; default 100)
%
%   Output:
%     I - numerical approximation of the integral
%
%   Example:
%     I = simpsons(@(x) exp(-x.^2), 0, 1, 100)

    if nargin < 4, n = 100; end
    if mod(n, 2) ~= 0
        n = n + 1;
        warning('simpsons: n must be even; using n = %d.', n);
    end

    h = (b - a) / n;
    x = a:h:b;
    y = arrayfun(f, x);

    I = h/3 * (y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end));
end
