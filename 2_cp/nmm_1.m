%% NMM

format long;

% % zadani z hodiny
% max_iter = 1e3;
% A = [7,6,-3;-12,-20,24;-6,-12,16];
% x0 = [1;2;3];
% tol = 1e-5;

% ukol 1.
% alpha = 30;  % a)
% alpha = -30;  % b)
% A = [alpha, 2, 3, 13; 5, 11, 10, 8; 9, 7, 6, 12; 4, 14, 15, 1];
% tol = 1e-6;
% x0 = [1;2;3;4];
% max_iter = 1e3;

% 4. ukol
A = [1, -2, -3, 0; 5, -9, 4, -3; 0, 0, -6, 1; -1, 3, 4, 6];
tol = 1e-5;
x0 = [1;2;3;4];
max_iter = 1e3;

function [x_out, sigma_out] = nmm(max_iter, A, x0, tol)
    eig_n_tol = tol;
    eig_v_tol = eig_n_tol;
    init_sigma = 1e10;
    
    % normalize x
    x0 = x0 / norm(x0);
    
    % global variables
    y_array = zeros(numel(x0), max_iter);
    x_array = zeros(numel(x0), max_iter);
    sigma_array = NaN(1, max_iter);
    x_array(:,1) = x0/norm(x0);
    sigma_array(1) = init_sigma;
    iter_cnt = 0;
    x_out = 0;
    sigma_out = 0;
    
    % indexing init + 1
    for i = 2:max_iter
        iter_cnt = i - 1;
        y_array(:,i) = A * x_array(:,i-1);
        x_array(:,i) = y_array(:,i) / norm(y_array(:,i));
        sigma_array(i) = x_array(:,i)' *  A * x_array(:,i);
    
        if abs(sigma_array(i)-sigma_array(i-1)) < eig_n_tol && ...
           norm(abs(x_array(:,i))-abs(x_array(:,i-1))) < eig_v_tol
            fprintf("Dosazeno pozadovane presnosti (iter_cnt=%d), breaking loop.\n",iter_cnt);
            x_out = x_array(:,i);
            sigma_out = sigma_array(i);
            break;
        end
    
        if i == max_iter
            fprintf("Dosazeno maxima iteraci, exiting loop.\n\n");
        end
    end
    
    % comparison
    [V, lambda] = eigs(A,1);
    disp("srovnani 1:"+abs(sigma_out - lambda)+"; sigma_out: "+sigma_out+ ...
         "; lambda: "+lambda);
    disp("srovnani 2:"+max(abs(V) - abs(x_out)));
end

[V, lambda] = nmm(max_iter, A, x0, tol)