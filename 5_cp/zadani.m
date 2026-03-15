%% Gaussova explicitni metoda
% Ukol a)
% func call
f = @(t,y) t*y+t^3;

[t,y] = cp5_expl_gauss(f, [0 1], 1, 10)

% Ukol b)
% call the function 
f = @(t,y)[y(2)^2 - 2*y(1); y(1) - y(2) - t*y(2)^2];
tspan = [0 1];
y0 = [0 1];
N = 10;

[t,y] = cp5_expl_gauss(f, tspan, y0, N)

%% Gaussova implicitni metoda 
% Ukol a)
% func call
f = @(t,y) t*y+t^3;

[t,y] = cp5_impl_gauss(f, [0 1], 1, 10)

% Ukol b)
% call the function 
f = @(t,y)[y(2)^2 - 2*y(1); y(1) - y(2) - t*y(2)^2];
tspan = [0 1];
y0 = [0 1];
N = 10;

[t,y] = cp5_impl_gauss(f, tspan, y0, N)

%% Gaussova lichobeznikova metoda
% NOTE:
% priklad hledani analytickeho reseni:
% f=@(t,y)(2*y+t)/t; tspan[1,2]; y0=1;
% syms y(t)
% presne = matlabFunction(dsolve(diff(y) == f(t,y), y(tspan(1)) == y0));
% clear y(t)

% Ukol a)
% func call
f = @(t,y) t*y+t^3;

[t,y] = cp5_lichb_gauss(f, [0 1], 1, 10)

% Ukol b)
% call the function 
f = @(t,y)[y(2)^2 - 2*y(1); y(1) - y(2) - t*y(2)^2];
tspan = [0 1];
y0 = [0 1];
N = 10;

[t,y] = cp5_lichb_gauss(f, tspan, y0, N)