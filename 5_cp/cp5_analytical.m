% function [t,yA,y_exact_sym]=cp5_analytical(f_sym,tspan,y0,N)
% 
%     t = linspace(tspan(1), tspan(2), N+1);
% 
%     % detect variables used in f_sym
%     vars = symvar(f_sym);
%     tt = vars(1);     % independent variable
%     yy = vars(2);     % state variable
% 
%     syms y(t_sym)
% 
%     % substitute variables properly
%     ode = diff(y,t_sym) == subs(f_sym, [tt, yy], [t_sym, y]);
% 
%     % solve ODE
%     y_exact_sym = dsolve(ode, y(tspan(1)) == y0);
% 
%     % numeric function
%     y_exact = matlabFunction(y_exact_sym, 'Vars', t_sym);
%     yA = y_exact(t);
% 
% end



%% if only y is present (no t) :)))) i love matlab
function [t,yA,y_exact_sym]=cp5_analytical(f_sym,tspan,y0,N)

    t = linspace(tspan(1), tspan(2), N+1);

    vars = symvar(f_sym);

    if length(vars) == 1
        yy = vars(1);      % only y present
        syms tt
    else
        tt = vars(1);
        yy = vars(2);
    end

    syms y(t_sym)

    ode = diff(y,t_sym) == subs(f_sym, [tt, yy], [t_sym, y]);

    y_exact_sym = dsolve(ode, y(tspan(1)) == y0);

    y_exact = matlabFunction(y_exact_sym, 'Vars', t_sym);
    yA = y_exact(t);

end