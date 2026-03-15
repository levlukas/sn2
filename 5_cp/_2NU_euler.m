%% zadani
x_0 = [0;0;0];

A_1 = [ 3 -1  1;
      1 -8 -2;
      1  1  5];
b_1 = [-2; 1 ; 4];


A_2 = [5 -7;
      -7 10];
b_2 = [-2 ; 3];

function D = d_mat(A)
    % A ... vstupni matice
    D = zeros(size(A));
    for i = 1:size(A,1)
        D(i,i) = A(i,i);
    end
end

function U = u_mat(A)
    % A ... vstupni matice
    U = zeros(size(A));
    for i = 1:size(A,1)
        for j = i+1:size(A,2)
            U(i,j) = A(i,j);
        end
    end
end

function L = l_mat(A)
    % A ... vstupni matice
    L = zeros(size(A));
    for i = 1:size(A,1)
        for j = 1:i-1
            L(i,j) = A(i,j);
        end
    end
end


function [x, iter_counter] = gs_iter(x_0, A, b, epsilon_threshold)
    % Resi soustavu linearnich rovnic Gauss-Seidelovou iteracni metodou.
    % x_0 ... pocatecni aproximace
    % A ... vstupni matice
    % b ... vstupni vektor pravych stran
    % epsilon_threshold ... kriterium pro zastaveni iterace
    iter_counter = 0;
    stop = false;  
    U = u_mat(A);
    D = d_mat(A);
    L = l_mat(A);
    LD = L+D;
    x_prev = x_0;
    try
        while not(stop)
            iter_counter = iter_counter + 1;
            x_latest = LD \ (b - U * x_prev);
            stop = norm(x_latest - x_prev,1) <= epsilon_threshold;  % pozije manhattan normu
            x_prev = x_latest;  % update pomocne promenne
        end
        x = x_latest;  % deklarace vystupu
    catch
        x = NaN;  % pokud nastala chyba, vrat nan
    end
end

%{
disp(d_mat(A_1));  % debug
disp(u_mat(A_1));  % debug
disp(l_mat(A_1));  % debug
disp(d_mat(A_1)+l_mat(A_1))  % debug
%}

% reseni soustavy 1
fprintf("==================\nReseni soustavy c. 1:\n");
[res_1, counter_1] = gs_iter(x_0, A_1, b_1, 1e-8);
disp(res_1);
fprintf("Pocet provedenych kroku: %d\n\n\n", counter_1);

% reseni soustavy 2 s ostre diag. nedominantni matici
fprintf("==================\nReseni soustavy c. 2:\n");
[res_2, counter_2] = gs_iter([0;0], A_2, b_2, 1e-3);
if isnan(res_2)
    fprintf("Chyba pri reseni matice A_2, pravdepodobne nekonverguje.\n" + ...
            "Pocet kroku pred dosazenim chyby: %d\n\n\n", counter_2);
else
    disp(res_2);
    fprintf("Pocet provedenych kroku: %d", counter_2);
end