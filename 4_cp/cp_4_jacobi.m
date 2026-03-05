function [VC, VV] = jacobi(A, tol, MaxIt)
%------------------------------------------------------------
% VSTUPY: 
% A=        % matice A
% tol =     % zadana presnost (nulovost mimodiag. prvku)
% MaxIt=    % nastaveni max. poctu iteraci
%------------------------------------------------------------
format long

% overte, zda A je symetricka, pokud neni ukoncete vypocet!
if issymmetric(A) == false
    return;
end
    
Ak=A;  
n = size(Ak,1);
Js = eye(n);   % jednotkova matice velikosti A
for k=1:MaxIt
    % najdete radek pk a sloupec qk v matici Ak (napr. pomoci prikazu max), 
    %   coz je pozice v absolutni hodnote nejvetsiho naddiagonalniho prvku 
    %   v matici Ak
    AU = triu(Ak,1);
    [~,ind] = max(abs(AU(:)));
    [pk,qk] = ind2sub(size(AU),ind);

    if abs(Ak(pk,qk))< tol, break, end % podminka ukonceni
        tau = (Ak(qk,qk)-Ak(pk,pk)) / (2 * Ak(qk,pk));  % spocteme tau
        if tau==0  % kvuli matlab fci sign, kde sign(0) = 0
            t=1;
        else
        t = sign(tau) / (abs(tau) + sqrt(1 + tau^2)); % spocteme t
    end
    c = 1 / (sqrt(1 + t^2)); % spocteme c
    s = c * t; % spocteme s

    % vypocet Jk lze udelat vypocetne jednoduseji bez vyjadreni cele ridke
    %   matice J
    Jk=eye(n); % jednotkova matice velikosti A 
    Jk([pk,qk],[pk,qk]) = [c s; -s c];  % zmente prvky v Jk na pozicich pp, qq, pq a qp na c a +-s
    Ak = Jk' * Ak * Jk;  % nova matice Ak pomoci matice Jk
    Js = Js * Jk;  % prinasobte k matici Js dalsi matici Jk 
end

if k==MaxIt, disp('dosazeno maxit!!!!'), end

%------------------------------------------------------------
% VYSTUPY:
[VC,idx] = sort(diag(Ak))  % vektor vlastnich cisel (jsou na diagonale Ak), serazeny prikazem sort
VV = Js(:,idx)  % matice vlastnich vektoru (vytvorily se v Js), serazeny dle poradi VC
PocetIteraci = k  % pocet iteraci  
%------------------------------------------------------------

%--- srovnání s EIG -----
[VVeig,VCeig]=eig(A); 
[VCeig,poradi]=sort(diag(VCeig));
VVeig2=VVeig(:,[poradi]);
fprintf('rozdil vl. cisel: %e\n',abs(VCeig-VC))
fprintf('maximalni rozdil vl. vektoru: %e\n',max(abs(abs(VVeig2)-abs(VV))))
end  % end function


A=[2.25 -0.25 -1.25 2.75;-0.25 2.25 2.75 1.25;-1.25 2.75 2.25 -0.25;2.75 1.25 -0.25 2.25];
tol=1e-5;
MaxIt= 1e5;

[VC, VV] = jacobi(A, tol, MaxIt);