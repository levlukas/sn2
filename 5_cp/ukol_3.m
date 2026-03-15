%% global variables
% N - vybirano tak, aby chyba dostatecne mala
N = 40;

g = 9.81;
l = 1;
f = @(t,y) [y(2); -g / l * sin(y(1))];

y0 = [0.5 * pi; 0];
tspan = [0 10];

for N = 1:60
    %% call
    
    [t_e, y_e] = cp5_expl_gauss(f, tspan, y0, N);
    [t_l, y_l] = cp5_lichb_gauss(f, tspan, y0, N);
    
    S = ode45(f, tspan, y0);
    
    fprintf("`ode45` provedlo %d kroku.\n", length(S.x));
    
    %% save figure
    folder = "ukol_3_plots";
    % create folder if it does not exist
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    filename = sprintf("N_%d.png", N);
    filepath = fullfile(folder, filename);
    saveas(gcf, filepath);
end