function [VC, VV] = bisection(A, VCk, VCkj, tol)
%            METODA BISEKCE
%------------------------------------------------------------
% % VSTUPY: 
% A =     % matice A
% VCk =   % od kolikateho vlastniho cisla chci hledat (index)
% VCkj =  % po kolikate vlastniho cislo chci hledat (index)
% tol =   % zadana presnost
%------------------------------------------------------------
format long

% ověření zda je matice A je symetrická (pokud neni koncim)
if issymmetric(A) == false
    return;
end

% zjisteni zda je matice A tridiagonalni, pokud neni tak upravit pomoci A=hess(A)
% reseni: pomoci vytazeni horniho trojuhelniku o spravnem rozmeru (uz je
%   symatricka A, takze staci overit)
if any(triu(A,3) ~= 0, 'all')  % nejde pouzit all - jiny prikaz
    A = hess(A);
end

aa = diag(A); % vytahnu z A jeji diagonalni prvky 
bb = diag(A,-1);  %  vytahnu z A jeji poddiagonalni prvky 
n = size(A,1);  % pocet radku matice A

% pozn: char. polynomy maji Sturmovu vlastnost (viz sesit)
% tvorba polynomu p:
p = cell(1,n+1);                
p{1} = @(x)1 ;                  % polynom p_0(x)
p{2} = @(x)(aa(1) - x);         % polynom p_1(x)

% vyplneni arraye polynomu
for i = 3:n+1                   % polynomy p_2(x),...,p_(n+1)(x)
    p{i} = @(x)(aa(i-1)-x).*p{i-1}(x) - (bb(i-2).^2)*p{i-2}(x);
end

% doplneni bb o 2 nuly kvuli indexaci
bb = vertcat(0, bb, 0);

% TODO: chyba
% najdete zacatek intervalu <a,b>, kde se nachazi vsechna vl. cisla matice A
a0 = min(aa(1:n) - (abs(bb(2:n+1)) + abs(bb(1:n))));
% najdete konec intervalu <a,b>, kde se nachazi vsechna vl. cisla matice A
b0 = max(aa + (abs(bb(2:n+1)) + abs(bb(1:n))));

VC = zeros(1,VCkj-VCk+1);  % deklarace vektoru VC

% Vypocet jednotlivych vl. cisel za pomoci Sturmovy vlastnosti
a=a0;
b=b0;
for i = 1:(VCkj-VCk+1)    
    while (b-a) > tol
      c = mean([a,b]);  % najdu stred intervalu (a,b) = (ak, bk)
      PocetZnamZmen = 0;
        for j=1:n
            if p{j}(c)*p{j+1}(c)<=0
                PocetZnamZmen = PocetZnamZmen+1;      % Pocet znamenkovych zmen v c
            end
        end
        if PocetZnamZmen < (VCkj-i+1)
            a = c;  % zmensim interval smerem doprava
        else
            b = c;  % zmensim interval smerem doleva
        end
    end
    b=0.5*(a+b);  % dopresnim nalezene vl. cislo
    a=a0;
    VC(i)=b;
end
    
%------------------------------------------------    
% VYSTUPY:
% VC = vektor hledanyvh vybranych vlastnich cisel    
VC=sort(VC);
%------------------------------------------------    

% Vykresleni polynomu p{n+1} a y=0
figure(1)
t=a0:0.01:b0;
plot(t,p{n+1}(t),'g');
hold on;
plot(t,0*t,'k');
hold on;
plot(VC,0*VC,'ro');

%--- srovnání s EIG -----
VCeig=eig(A); 
VCeig=sort(VCeig);
fprintf('rozdil vl. cisel: %e\n',abs(VCeig(VCk:VCkj) -VC'))

end  % end function

A=[2.25 -0.25 -1.25 2.75;-0.25 2.25 2.75 1.25;-1.25 2.75 2.25 -0.25;2.75 1.25 -0.25 2.25];
tol=10^(-5); 
VCk = 2;
VCkj = 3;

[VC, VV] = bisection(A, VCk, VCkj, tol)