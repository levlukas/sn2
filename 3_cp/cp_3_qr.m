
format long

% inputs
% 1) sym. matice
% A = [2.25 -0.25 -1.25 2.75;-0.25 2.25 2.75 1.25;-1.25 2.75 2.25 -0.25;...
     % 2.75 1.25 -0.25 2.25];
% 2) nesym. matice
A = [2.5 -2.5 3 0.5;0 5 -2 2;-0.5 -0.5 4 2.5;-2.5 -2.5 5 3.5];

% Hessenberguv tvar je treba pro nektere z uloh (1b), 2b))
A = hess(A);

%% VSTUPY: 
% - matice A
% - tol = zadana tolerance
% - MaxIter = maximalni pocet iteraci
% VYSTUPY:
% VC = vektor vlastnich cisel
% VV = matice vlastnich vektoru
% k = pocet iteraci
%%------------------------------------------------------------
function [v1] = INMMsPgrader(A,mu,x0,tol,MaxIter)
    max_iter = MaxIter;
    inp_mat = A;
    inp_vec = x0;
    eig_n_tol = tol;
    eig_v_tol = tol;
    init_sigma = 1e10;
    % normalize x
    inp_vec = inp_vec / norm(inp_vec);
    % global variables
    y_array = zeros(numel(inp_vec), max_iter);
    x_array = zeros(numel(inp_vec), max_iter);
    sigma_array = NaN(1, max_iter);
    x_array(:,1) = inp_vec/norm(inp_vec);
    sigma_array(1) = init_sigma;
    iter_cnt = 0;
    x_out = 0;
    sigma_out = 0;
    % special variable (method exclusive)
    [L,U,P] = lu(inp_mat - mu * eye(size(inp_mat)));
    % indexing init + 1
    for i = 2:max_iter
        iter_cnt = i - 1;
        z_temp = L \ (P * x_array(:,i-1));
        y_array(:,i) = U \ z_temp;
        x_array(:,i) = y_array(:,i) / norm(y_array(:,i));
        sigma_array(i) = x_array(:,i)' *  inp_mat * x_array(:,i);
        if norm(abs(x_array(:,i))-abs(x_array(:,i-1))) < eig_v_tol
            fprintf("Dosazeno pozadovane presnosti (iter_cnt=%d), breaking loop.\n",iter_cnt);
            x_out = x_array(:,i);
            sigma_out = sigma_array(i);
            break;
        end
        if i == max_iter
            fprintf("Dosazeno maxima iteraci, exiting loop.\n\n");
        end
    end
    v1 = x_out;
end

function [VC, VV, k]=QRalgoritmus(A,tol,MaxIter)
    % assert: A = trans(A)
    if A == A'
        symetr=1;
    else
        symetr=0;
    end
    
    % init Q (cummulative)
    Q = eye;
    Ak=A;  % A vlozim do Ak
    for k=1:MaxIter
        [Qk,Rk] = qr(Ak);  % QR rozklad z A_k-1
         Ak= Rk * Qk;  % vypocet A_k
         Q = Q * Qk;  % vytvarim postupne matici Q=Q1*Q2*Q3...*Qk  
            if abs(tril(Ak,-1)) < tol
                break;
            end
    end
    % hledana vlastni cisla jsou na diagonale Ak
    VC = diag(Ak);

    if k==MaxIter, fprintf("Dosazeno max. iteraci (%d)", MaxIter), end
    
    if symetr == true
         disp('matice A je symetricka'); 
         VV=Q;
    else       
        disp('matice A neni symetricka'); 
        for i=1:length(VC)
             [v1] = INMMsPgrader(A,VC(i)-.01,1:length(VC),tol,MaxIter);
             VV(:,i)=v1; % zapiseme nalezeny vektor do matice vlastnich vektoru VV
        end
    end
    
    % func call
    [VC, VV, k] = QRalgoritmus(A,1e-6,1e5)
    
    %--- srovnání s EIG -----
    [VVeig,VCeig]=eig(A); 
    [VCeig,poradi]=sort(diag(VCeig));
    [VCmoje,poradi2]=sort(VC);
    VVeig2=VVeig(:,[poradi]);
    VVmoje=VV(:,[poradi2]);
    fprintf('rozdil vl. cisel: %e\n',abs(VCeig-VCmoje))
    fprintf('maximalni rozdil vl. vektoru: %e\n',max(abs(abs(VVeig2)-abs(VVmoje))))
end  % func. def end