%% global variables
f = @(t,y) -13*y;  % rovnice
y0 = 1;  % poc. podminka
tspan = [0 1];
syms t_a y_a;
f_a = -13*y_a + 0*t_a;

outdir = "plots";
if ~exist(outdir,'dir')
    mkdir(outdir);
end

for N = [30 20 10 5]

    figure('Visible','off','Name',sprintf('N=%d',N),'NumberTitle','off');
    hold on; grid on;

    [t, y] = cp5_expl_gauss(f, tspan, y0, N);
    plot(t, y(1,:), '-o', 'DisplayName', 'Explicit Euler');

    [t, y] = cp5_impl_gauss(f, tspan, y0, N);
    plot(t, y(1,:), '-s', 'DisplayName', 'Implicit Euler');

    [t, y] = cp5_lichb_gauss(f, tspan, y0, N);
    plot(t, y(1,:), '-^', 'DisplayName', 'Lichobeznik');

    [t, yA, y_symbolic] = cp5_analytical(f_a, tspan, y0, N);
    plot(t, yA, '-k', 'LineWidth', 1.5, 'DisplayName', 'Analytical');

    xlabel('t');
    ylabel('y(t)');
    title(sprintf('Solutions for N = %d', N));
    legend('Location','best');

    filename = sprintf('%s/solution_N%d.png', outdir, N);
    saveas(gcf, filename);

    close(gcf)
end