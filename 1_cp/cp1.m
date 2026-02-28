%% Cp1
% Note:
% - popis chovani vl. c. a v. pro ruzne transf. matic
%   je popsan v prednasce /1t/
%% Ukol 1

A=[-14, -30, 42; 24, 49, -66; 12, 24, -32];

%% 1.1
[V, lambda] = eig(A)

%% 1.2
B = A - 3*eye;
[V, lambda] = eig(B)
% vektory jsou porad stejne, lambdy ne 

%% 1.3
C = inv(A);
[V, lambda] = eig(C)
% lambda_a^(-1) = lambda_c
% V_a = V_c

%% 1.4
D = A^(3);
[V, lambda] = eig(D)
% lambda_a^(3) = lambda_d
% V_a = V_d

%% 1.5
T=[1, 2, -1; 2, 0, 3;-1, 2, 2];
F = inv(T)*A*T;
[V, lambda] = eig(F)
% lambda_a = lambda_f
% V_f = transfromace z V_a

%% Ukol 2

%% 2.1
% sym. matrix
A=[5, 1, 1, 1, 1;
   1, 5, 1, 1, 1;
   1, 1, 5, 1, 1;
   1, 1, 1, 5, 1;
   1, 1, 1, 1, 5];

[V, lambda] = eig(A);
% res ... vysledky skalovany faktorem 1e-15, coz je ocekavany 'sum'
res = V'*V - eye(size(V))
% vektory jsou ortonormalni, takze velikost je 1
% ortonormalita rika, ze Q^T*Q = I (P/2t/)

%% 2.2
A=[-2, 5, 1, -3; 5, 4, -1, 2; 1, -1, 2, 0; -3, 2, 0, -3];
[V, lambda] = eig(A);
% nefunguje -> isreal(inp::matrix)
%for i=lambda(:)
%    if not isreal(i)
%        fprintf('Eigenvalue %.2f is not real\n', i);
%    end;
%end;
isreal(V)

%% 2.3
A = [1, -1, -1; 1, 1, 0; 3, 0, 1]
% A=[-7, 9; -1, -1]

[V, lambda_mat] = eig(A);
lambda = diag(lambda_mat);

% expected accuracy of lambda
lambda_tol = 1e-10;
lambda_rounded = round(lambda/lambda_tol)*lambda_tol;
unique_lambdas = unique(lambda_rounded);

% sanity check - defectivity
is_defective = false;
for k = 1:length(unique_lambdas)
    lam = unique_lambdas(k);

    % algebraic multiplicity
    alg_mult = sum(abs(lambda_rounded - lam) < lambda_tol);

    % geometric multiplicity
    geom_mult = size(A,1) - rank(A - lam*eye(size(A)));

    if geom_mult < alg_mult
        is_defective = true;
        break;
    end
end

if is_defective
    fprintf("Defektni matice! Abort\n");
    return;
end

D = inv(V)*A*V

%% 3
A=[1,-2,-3,0; 5,-9,4,-3; 0,0,-6,1; -1,3,4,6]

[V, lambda_mat] = eig(A);

lambda = diag(lambda_mat);

A_diagless = A-diag(A);
radii = zeros(size(A,1),1);
centers = zeros(size(A,1),1);
ax = axes;
for i = 1:size(A,1)
    centers(i) = A(i,i);
    radii(i) = sum(abs(A(i,:))) - abs(A(i,i));
end

scatter(ax, real(lambda), imag(lambda))
viscircles(ax,[real(centers), imag(centers)],radii)