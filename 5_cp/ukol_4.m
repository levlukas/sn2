%% global variables
f = @(t,y) t - t.*y/2;
y0 = 1;
tspan = [1 4];
N_arr = [100, 200, 300, 400];
tau_arr = (tspan(2) - tspan(1)) ./ N_arr;

syms t_sym y_sym
f_sym = t_sym - t_sym*y_sym/2;

E_expl = zeros(size(N_arr));
E_impl = zeros(size(N_arr));
E_lichb = zeros(size(N_arr));

%% call
for i = 1:numel(N_arr)
    N = N_arr(i);

    [~, yA, y_exact_sym] = cp5_analytical(f_sym, tspan, y0, N);
    [~, y_expl] = cp5_expl_gauss(f, tspan, y0, N);
    [~, y_impl] = cp5_impl_gauss(f, tspan, y0, N);
    [~, y_lichb] = cp5_lichb_gauss(f, tspan, y0, N);

    % globalni chyba v uzlech: E_N = max_n |y(t_n) - y_n|
    E_expl(i) = max(abs(y_expl(1,:) - yA));
    E_impl(i) = max(abs(y_impl(1,:) - yA));
    E_lichb(i) = max(abs(y_lichb(1,:) - yA));
end

% rozdily chyb mezi po sobe jdoucimi N (informativne)
dE_expl = E_expl(1:end-1) - E_expl(2:end);
dE_impl = E_impl(1:end-1) - E_impl(2:end);
dE_lichb = E_lichb(1:end-1) - E_lichb(2:end);

% odhad radu konvergence:
% E_N ~ C * tau^p  =>  p ~ log(E_N / E_M) / log(tau_N / tau_M)
p_expl = log(E_expl(1:end-1) ./ E_expl(2:end)) ./ ...
         log(tau_arr(1:end-1) ./ tau_arr(2:end));
p_impl = log(E_impl(1:end-1) ./ E_impl(2:end)) ./ ...
         log(tau_arr(1:end-1) ./ tau_arr(2:end));
p_lichb = log(E_lichb(1:end-1) ./ E_lichb(2:end)) ./ ...
          log(tau_arr(1:end-1) ./ tau_arr(2:end));

fprintf('Presne reseni: y(t) = %s\n\n', char(y_exact_sym));

fprintf('Chyby E_N = max_n |y(t_n) - y_n|\n');
fprintf('%6s %12s %16s %16s %16s\n', 'N', 'tau', 'E_explicit', 'E_implicit', 'E_lichb');
for i = 1:numel(N_arr)
    fprintf('%6d %12.6f %16.6e %16.6e %16.6e\n', ...
    N_arr(i), tau_arr(i), E_expl(i), E_impl(i), E_lichb(i));
end

fprintf('\nRozdily E_N - E_M pro po sobe jdoucich N (jen informativne)\n');
fprintf('%12s %16s %16s %16s\n', 'par (N,M)', 'dE_explicit', 'dE_implicit', 'dE_lichb');
for i = 1:numel(N_arr)-1
    fprintf('(%3d,%3d) %16.6e %16.6e %16.6e\n', ...
        N_arr(i), N_arr(i+1), dE_expl(i), dE_impl(i), dE_lichb(i));
end

fprintf('\nOdhad radu konvergence p\n');
fprintf('%12s %16s %16s %16s\n', 'par (N,M)', 'p_explicit', 'p_implicit', 'p_lichb');
for i = 1:numel(N_arr)-1
    fprintf('(%3d,%3d) %16.6f %16.6f %16.6f\n', ...
        N_arr(i), N_arr(i+1), p_expl(i), p_impl(i), p_lichb(i));
end