function [L, U, P] = lu_decomposition(A)
% LU_DECOMPOSITION  Compute the LU decomposition of A with partial pivoting
%                   such that P*A = L*U.
%
%   [L, U, P] = lu_decomposition(A)
%
%   Inputs:
%     A - n×n square matrix
%
%   Outputs:
%     L - lower triangular matrix (unit diagonal)
%     U - upper triangular matrix
%     P - permutation matrix
%
%   To solve Ax = b after decomposition:
%     y = L \ (P*b);
%     x = U \ y;
%
%   Example:
%     A = [2 1 -1; -3 -1 2; -2 1 2];
%     [L, U, P] = lu_decomposition(A)
%     b = [8; -11; -3];
%     x = U \ (L \ (P*b))   % should give [2; 3; -1]

    n = size(A, 1);
    U = A;
    L = eye(n);
    P = eye(n);

    for k = 1:n-1
        % Partial pivoting
        [~, p] = max(abs(U(k:n, k)));
        p = p + k - 1;

        if p ~= k
            U([k, p], :)     = U([p, k], :);
            P([k, p], :)     = P([p, k], :);
            L([k, p], 1:k-1) = L([p, k], 1:k-1);
        end

        if abs(U(k,k)) < eps
            error('lu_decomposition: matrix is singular or nearly singular.');
        end

        for i = k+1:n
            L(i,k) = U(i,k) / U(k,k);
            U(i,k:n) = U(i,k:n) - L(i,k) * U(k,k:n);
        end
    end
end
