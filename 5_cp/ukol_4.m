%% global variables
f = @(t,y) t - t*y/2;
y0 = 1;
tspan = [1 4];
N_arr = [100, 200, 300, 400];
tau_arr = N_arr ./ (tspan(2) - tspan(1));
E_N = zeros(1,size(N_arr,2)-1);

%% call
for i = 1:4
    cp5_expl_gauss(f, tspan, y0, N_arr(i));
    cp5_impl_gauss(f, tspan, y0, N_arr(i));
    cp5_lichb_gauss(f, tspan, y0, N_arr(i));

    if i > 1
        
    end
end