function [t,y]=cp5_impl_gauss(f,tspan,y0,N)
    % implicitní Eulerova metoda

    % rovnomerne rozdeleni intervalu tspan=<t0,tN> na N+1 uzlovych bodu t0,t1,..,tN
    t = linspace(tspan(1), tspan(2), N+1);
    % delka casoveho kroku
    tau = t(2) - t(1);  % t je ekvidist.

    y0 = y0(:);          % vzdy sloupcove
    y = zeros(numel(y0), N+1);  % predalokace matice reseni

    % nacteni pocatecni podminky do matice reseni
    y(:,1) = y0;

    for n=1:N
        % vypocet reseni v novem casovem kroku, jako řešič nelineární
        % soustavy použijte fsolve
        % pouziti fsolve: `fsolve(@(z) F, z0)`, kde z je obecne nejaka random
        % variable, F je soustava zpasana ve tvaru F=0, poc. iterace je z0

        % y(:,n+1) - y(:,n) - tau * f(t(n+1),y(:,n+1)) = F = 0
        % moje neznama je ale vektor, takze substituce y(:,n+1) = z
        % z - y(:,n) - tau * f(t(n+1),z) = F = 0

        % protoze resim pro stupen n+1 (indexace mtlb od 1), pak lze zapsat
        % jednoradkove:
        % y(:,n+1) = fsolve(@(z) F(z), z0)

        F  = @(z) z - y(:,n) - tau * f(t(n+1), z);  % soustava rovnic F(z)=0
        z0 = y(:,n);                                % rozumny pocatecni odhad (predchozi krok)

        % volitelne: potlacit vypis z fsolve
        opts = optimoptions('fsolve','Display','off');

        y(:,n+1) = fsolve(F, z0, opts);
    end
end

