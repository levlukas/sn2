function I = trapezoidal(f, a, b, n)
% TRAPEZOIDAL  Approximate the integral of f over [a,b] using the
%              composite trapezoidal rule.
%
%   I = trapezoidal(f, a, b, n)
%
%   Inputs:
%     f    - function handle
%     a, b - integration limits
%     n    - number of sub-intervals (default 100)
%
%   Output:
%     I - numerical approximation of the integral
%
%   Example:
%     I = trapezoidal(@(x) sin(x), 0, pi, 1000)
%     % returns approximately 2.0000

    if nargin < 4, n = 100; end

    h = (b - a) / n;
    x = a:h:b;
    y = arrayfun(f, x);

    I = h * (y(1)/2 + sum(y(2:end-1)) + y(end)/2);
end
