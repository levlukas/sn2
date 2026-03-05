function [df_fwd, df_bwd, df_cen] = finite_diff(f, x, h)
% FINITE_DIFF  Estimate the first derivative of f at x using finite
%              difference formulas.
%
%   [df_fwd, df_bwd, df_cen] = finite_diff(f, x, h)
%
%   Inputs:
%     f - function handle
%     x - point at which to differentiate
%     h - step size (default 1e-5)
%
%   Outputs:
%     df_fwd - forward  difference: [f(x+h) - f(x)  ] / h        O(h)
%     df_bwd - backward difference: [f(x)   - f(x-h)] / h        O(h)
%     df_cen - central  difference: [f(x+h) - f(x-h)] / (2h)    O(h^2)
%
%   Example:
%     [fwd, bwd, cen] = finite_diff(@(x) sin(x), pi/4, 1e-5)
%     % all three should be close to cos(pi/4) ≈ 0.7071

    if nargin < 3, h = 1e-5; end

    df_fwd = (f(x + h) - f(x)    ) / h;
    df_bwd = (f(x)     - f(x - h)) / h;
    df_cen = (f(x + h) - f(x - h)) / (2 * h);
end
