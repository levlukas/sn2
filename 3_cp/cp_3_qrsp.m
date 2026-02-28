%% init values
% 1. ukol
% A=[2.25 -0.25 -1.25 2.75;-0.25 2.25 2.75 1.25;-1.25 2.75 2.25 -0.25;2.75 1.25 -0.25 2.25];
% 1b.
% A = hess(A);

% % 2.
% A = [2.5 -2.5 3 0.5;0 5 -2 2;-0.5 -0.5 4 2.5;-2.5 -2.5 5 3.5];
% % 2b.
% A = hess(A);

tol = 1e-6;
MaxIter = 1e5;

function [VC,kk]=QRalgoritmusPosun(A,tol,MaxIter)
    %% VSTUPY: 
    % - matice A
    % - tol = zadana tolerance
    % - MaxIter = maximalni pocet iteraci
    % VYSTUPY:
    % VC = vektor vlastnich cisel
    % k = pocet iteraci
    %%------------------------------------------------------------
    
    format long
    
    Ak=A;  % A vlozim do Ak
    
    n=size(Ak,1);
    
    kk=0;  % celkovy pocet iteraci
    for k=n:-1:2
        I=eye(k);
        Ak=Ak(1:k,1:k);  % z matice Ak vyberu submatici s kxk
        for i=1:MaxIter
            muk=Ak(k,k);  % zvolime Rayleighuv posun
            % provedeme QR rozklad posunute matice Ak o muk (r.4 v alg. 7.8)
            [Qk,Rk] = qr(Ak - muk*I);  
            Ak = Rk * Qk + muk * I;  % zmenime Ak dle r.5 v alg. 7.8
            kk=kk+1;  % zvysime citac celkovych iteraci
            
            % podm. ukonceni: |vsech prvku v poslednim radku matice Ak, mimo prvek A(k,k)|< tol  
            if all(abs(Ak(k,1:k-1)) < tol)
                break;
            end
        end
        VC(k,1)=Ak(k,k); % Ak(k,k) dokonvergovalo k vl. cislu 
    end
    VC(1,1)=Ak(1,1);
    
    if i==MaxIter disp('dosazeno MaxIt-> NEKONVERGUJE!'), return; end
    
    %--- srovnání s EIG -----
    [VCeig]=eig(A); 
    VCeig=sort(VCeig);
    VC=sort(VC);
    fprintf('rozdíl vl. cisel: %e\n',abs(VCeig-VC));
end

% function call
[VC,kk]=QRalgoritmusPosun(A,tol,MaxIter)