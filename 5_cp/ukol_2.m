%% glob variables
tspan = [1 3];
f = @(y,t) t*y/(t-2); 
y0 = 1;
syms y_a t_a
f_a = t_a * y_a / (t_a - 2);
N = 25;

%% call
S_23  = ode23(f, tspan, y0);
S_78  = ode78(f, tspan, y0);
S_113 = ode113(f, tspan, y0);
[t,y] = cp5_expl_gauss(f,tspan,y0,N);

%% plot all varibales

% create folder if it does not exist
folder = "ukol_2_plots";
if ~exist(folder, 'dir')
    mkdir(folder);
end

figure;
hold on;
plot(S_23.x, S_23.y, 'r-+', 'DisplayName', 'ode23');
plot(S_78.x, S_78.y, 'g-+', 'DisplayName', 'ode78');
plot(S_113.x, S_113.y, 'b-+', 'DisplayName', 'ode113');
plot(t, y, 'k-+', 'DisplayName', 'cp5\_expl\_gauss');
xlabel('Time');
ylabel('Solution y');
legend show;
hold off;

filename = sprintf("ukol_2.png");
filepath = fullfile(folder, filename);

saveas(gcf, filepath);