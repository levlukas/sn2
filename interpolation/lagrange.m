function yq = lagrange(x, y, xq)
% LAGRANGE  Evaluate the Lagrange interpolating polynomial at points xq.
%
%   yq = lagrange(x, y, xq)
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
%     x = [0, 1, 2];
%     y = [1, 3, 2];
%     yq = lagrange(x, y, 1.5)   % should return 2.875

    n  = length(x);
    nq = length(xq);
    yq = zeros(1, nq);

    for q = 1:nq
        s = 0;
        for i = 1:n
            % Build the i-th Lagrange basis polynomial evaluated at xq(q)
            Li = 1;
            for j = 1:n
                if j ~= i
                    Li = Li * (xq(q) - x(j)) / (x(i) - x(j));
                end
            end
            s = s + y(i) * Li;
        end
        yq(q) = s;
    end
end
