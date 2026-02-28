function x = gauss_elimination(A, b)
% GAUSS_ELIMINATION  Solve the linear system Ax = b using Gaussian
%                    elimination with partial (row) pivoting.
%
%   x = gauss_elimination(A, b)
%
%   Inputs:
%     A - n×n coefficient matrix
%     b - n×1 right-hand side vector
%
%   Output:
%     x - n×1 solution vector
%
%   Example:
%     A = [2 1 -1; -3 -1 2; -2 1 2];
%     b = [8; -11; -3];
%     x = gauss_elimination(A, b)
%     % returns [2; 3; -1]

    [n, ~] = size(A);
    Ab = [A, b];  % augmented matrix

    % Forward elimination with partial pivoting
    for k = 1:n-1
        % Pivot: find row with maximum absolute value in column k
        [~, p] = max(abs(Ab(k:n, k)));
        p = p + k - 1;
        if p ~= k
            Ab([k, p], :) = Ab([p, k], :);
        end

        if abs(Ab(k,k)) < eps
            error('gauss_elimination: matrix is singular or nearly singular.');
        end

        for i = k+1:n
            factor = Ab(i,k) / Ab(k,k);
            Ab(i, k:end) = Ab(i, k:end) - factor * Ab(k, k:end);
        end
    end

    % Back substitution
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (Ab(i, end) - Ab(i, i+1:n) * x(i+1:n)) / Ab(i,i);
    end
end
