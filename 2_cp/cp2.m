%% NMM

max_iter = 1e3;
inp_mat = [7,6,-3;-12,-20,24;-6,-12,16];
inp_vec = [1;2;3];
eig_n_tol = 1e-5;
eig_v_tol = eig_n_tol;
init_sigma = 1e10;

% normalize x
inp_vec = inp_vec / norm(inp_vec);

% global variables
y_array = zeros(numel(inp_vec), max_iter);
x_array = zeros(numel(inp_vec), max_iter);
sigma_array = NaN(1, max_iter);
x_array(:,1) = inp_vec/norm(inp_vec);
sigma_array(1) = init_sigma;
iter_cnt = 0;
x_out = 0;
sigma_out = 0;

% indexing init + 1
for i = 2:max_iter
    iter_cnt = i - 1;
    y_array(:,i) = inp_mat * x_array(:,i-1);
    x_array(:,i) = y_array(:,i) / norm(y_array(:,i));
    sigma_array(i) = x_array(:,i)' *  inp_mat * x_array(:,i);

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
[V, lambda] = eigs(inp_mat,1);
disp("srovnani 1:"+abs(sigma_out - lambda)+"; sigma_out: "+sigma_out+ ...
     "; lambda: "+lambda);
disp("srovnani 2:"+max(abs(V) - abs(x_out)));

%% INMM

max_iter = 1e3;
inp_mat = [7,6,-3;-12,-20,24;-6,-12,16];
inp_vec = [1;2;3];
eig_n_tol = 1e-5;
eig_v_tol = eig_n_tol;
init_sigma = 1e10;

% normalize x
inp_vec = inp_vec / norm(inp_vec);

% global variables
y_array = zeros(numel(inp_vec), max_iter);
x_array = zeros(numel(inp_vec), max_iter);
sigma_array = NaN(1, max_iter);
x_array(:,1) = inp_vec/norm(inp_vec);
sigma_array(1) = init_sigma;
iter_cnt = 0;
x_out = 0;
sigma_out = 0;

% special variable (method exclusive)
[L,U,P] = lu(inp_mat);

% indexing init + 1
for i = 2:max_iter
    iter_cnt = i - 1;
    z_temp = L \ (P * x_array(:,i-1));
    y_array(:,i) = U \ z_temp;
    x_array(:,i) = y_array(:,i) / norm(y_array(:,i));
    sigma_array(i) = x_array(:,i)' *  inp_mat * x_array(:,i);

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
[V, lambda] = eigs(inp_mat,1,'smallestabs');
disp("srovnani 1:"+abs(sigma_out - lambda)+"; sigma_out: "+sigma_out+ ...
     "; lambda: "+lambda);
disp("srovnani 2:"+max(abs(V) - abs(x_out)));

%% INMMsP

max_iter = 1e3;
inp_mat = [7,6,-3;-12,-20,24;-6,-12,16];
inp_vec = [1;2;3];
eig_n_tol = 1e-5;
eig_v_tol = eig_n_tol;
init_sigma = 1e10;
mu = -2.01;

% normalize x
inp_vec = inp_vec / norm(inp_vec);

% global variables
y_array = zeros(numel(inp_vec), max_iter);
x_array = zeros(numel(inp_vec), max_iter);
sigma_array = NaN(1, max_iter);
x_array(:,1) = inp_vec/norm(inp_vec);
sigma_array(1) = init_sigma;
iter_cnt = 0;
x_out = 0;
sigma_out = 0;

% special variable (method exclusive)
[L,U,P] = lu(inp_mat - mu * eye(size(inp_mat)));

% indexing init + 1
for i = 2:max_iter
    iter_cnt = i - 1;
    z_temp = L \ (P * x_array(:,i-1));
    y_array(:,i) = U \ z_temp;
    x_array(:,i) = y_array(:,i) / norm(y_array(:,i));
    sigma_array(i) = x_array(:,i)' *  inp_mat * x_array(:,i);

    if norm(abs(x_array(:,i))-abs(x_array(:,i-1))) < eig_v_tol
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
[V, lambda] = eigs(inp_mat, 1, mu);
disp("srovnani: "+max(abs(V) - abs(x_out)));