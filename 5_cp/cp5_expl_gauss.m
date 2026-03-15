function [t,y]=cp5_expl_gauss(f,tspan,y0,N)
    % explicitní Eulerova metoda

    % rovnomerne rozdeleni intervalu tspan=<t0,tN> na N+1 uzlovych bodu t0,t1,..,tN
    t = linspace(tspan(1), tspan(2), N+1);  
    tau = t(2) - t(1);  % delka casoveho kroku
    y0 = y0(:);  % vzdy sloupcove
    y(:,1)= y0; % nacteni pocatecni podminky do matice reseni

    for n=1:N
       % vypocet reseni v novem casovem kroku
       y(:,n+1) = y(:,n) + tau * f(t(n),y(:,n));
    end
end
