function [t,y]=cp5_lichb_gauss(f,tspan,y0,N)
    % cp5_lichb_gaussova metoda (trapezoidal rule / Crank–Nicolson)
    %
    % y_{n+1} = y_n + (tau/2) * ( f(t_n,y_n) + f(t_{n+1},y_{n+1}) )
    % => resime nelinearni soustavu pro y_{n+1} pomoci fsolve
    
    % rovnomerne rozdeleni intervalu tspan=<t0,tN> na N+1 uzlovych bodu t0,t1,..,tN
    t = linspace(tspan(1), tspan(2), N+1);  
    % delka casoveho kroku
    tau = t(2) - t(1);
    
    y0 = y0(:);          % vzdy sloupcove
    y = zeros(numel(y0), N+1);  % predalokace matice reseni

    % nacteni pocatecni podminky do matice reseni
    y(:,1) = y0;
    
    opts = optimoptions('fsolve','Display','off');
    
    for n=1:N
             % vypocet reseni v novem casovem kroku, jako řešič nelineární
             % soustavy použijte fsolve
             % F(z)=0, kde z = y(:,n+1):
             % z - y(:,n) - (tau/2)*( f(t(n),y(:,n)) + f(t(n+1),z) ) = 0
    
             Fn = f(t(n), y(:,n));
             F  = @(z) z - y(:,n) - (tau/2) * ( Fn + f(t(n+1), z) );
    
             % pocatecni odhad (predchozi krok)
             z0 = y(:,n);
             y(:,n+1) = fsolve(F, z0, opts);
    end
end