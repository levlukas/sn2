%% global variables
% N - vybirano tak, aby chyba dostatecne mala
N = 40;

g = 9.81;
l = 1;
f = @(t,y) [y(2); -g / l * sin(y(1))];

y0 = [0.5 * pi; 0];
tspan = [0 10];

folder = "ukol_3_plots";
if ~exist(folder, 'dir')
    mkdir(folder);
end

for N = 1:60
    [t_e, y_e] = cp5_expl_gauss(f, tspan, y0, N);
    [t_l, y_l] = cp5_lichb_gauss(f, tspan, y0, N);
    S = ode45(f, tspan, y0);

    % create new figure for this iteration (or clear existing)
    fig = figure('Visible','off'); % use 'on' to see during run
    hold on

    % plot only the first state (angle)
    plot(S.x,  S.y(1,:), 'b-+', 'DisplayName', 'ode45');
    plot(t_e,  y_e(1,:), 'r-+', 'DisplayName', 'expl');
    plot(t_l,  y_l(1,:), 'g-+', 'DisplayName', 'lichb');

    xlabel('Time'); ylabel('theta');
    legend('Location','best');
    hold off

    drawnow; % ensure graphics are rendered before saving

    filename = sprintf("N_%d.png", N);
    filepath = fullfile(folder, filename);
    saveas(fig, filepath);

    close(fig); % close figure to avoid too many open windows
end
