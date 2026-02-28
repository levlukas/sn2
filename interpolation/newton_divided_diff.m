function yq = newton_divided_diff(x, y, xq)
% NEWTON_DIVIDED_DIFF  Evaluate Newton's divided-difference interpolating
%                      polynomial at points xq.
%
%   yq = newton_divided_diff(x, y, xq)
%
%   Inputs:
%     x  - vector of n distinct interpolation nodes
%     y  - vector of n function values y(i) = f(x(i))
%     xq - scalar or vector of query points
%
%   Output:
%     yq - interpolated values at xq
%
%   Example:
%     x = [1, 3, 5, 7];
%     y = [3, 5, 4, 2];
%     yq = newton_divided_diff(x, y, [2, 4, 6])

    n = length(x);

    % Build divided-difference table; coefficients are the diagonal entries
    dd = y(:)';  % row vector copy
    coeff = zeros(1, n);
    coeff(1) = dd(1);

    for k = 2:n
        dd(k:n) = (dd(k:n) - dd(k-1:n-1)) ./ (x(k:n) - x(1:n-k+1));
        coeff(k) = dd(k);
    end

    % Evaluate using Horner-like nested multiplication
    nq = length(xq);
    yq = zeros(1, nq);
    for q = 1:nq
        val = coeff(n);
        for k = n-1:-1:1
            val = val * (xq(q) - x(k)) + coeff(k);
        end
        yq(q) = val;
    end
end
