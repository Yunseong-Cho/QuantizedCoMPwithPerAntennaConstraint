function [ lambda, count , SINR, D ] = Algo_joint_WB_PA(SP, G, n, initlambda, noPA)
% actual


b = SP.b;
if b == inf
    alpha = 1;
else
    bTable = [0.3634, 0.1175, 0.03454, 0.009497, 0.002499];
    if b > 5
        alpha = (1 - pi*sqrt(3)/2*2.^(-2*b));
    else
        alpha = (1 - bTable(floor(b)));
    end
end

gamma = SP.gamma(n);

Nr = SP.Nr;
Nc = SP.Nc;
Nu = SP.Nu;
th = SP.th;
Nsc = SP.Nsc;
algoMax = SP.algoMax;
lambda_curr = initlambda; % initial guess
lambda_past = zeros(Nu,Nc,Nsc);

W_DFT = dftmtx(Nsc)/sqrt(Nsc);
Psi = kron(W_DFT, eye(Nr));


% A = rand(Nr, Nr, Nc);

for i = 1:Nc
      D(:,:,i) = eye(Nr);
end
D_past = zeros(Nr, Nr, Nc);

count = 0;

for k = 1:Nsc
    Psi_k(:,:,k) = kron(W_DFT(k,:), eye(Nr));
end
% count2 = 1;
while(sum( abs((D(:) - D_past(:))./D(:)) > SP.th2) ~= 0)
    D_past = D;
    cnt = 0;
    
    count = count + 1;

    while (sum( abs((lambda_curr(:) - lambda_past(:))./lambda_curr(:)) > th) ~= 0)   % relative error, was abs(lambda_c - lambda_p)
    %while (1)
    cnt = cnt + 1;
    lambda_past = lambda_curr;
    
    
    LAMBDA = lambda_stack(lambda_past, SP);
    
    
        for k = 1:Nsc

            G_k = G{k}; % channel of k-th subcarrier 
            lambda_k = lambda_past(:,:,k);

            Lambda_k = diag(lambda_k(:));
            

            for i = 1:Nc
                G_i = G_stack(G,i,SP);

                ICI = G_k(:,:,i)*Lambda_k*G_k(:,:,i)';

                for u = 1:Nu                
                    Kz = alpha*ICI ...
                        + D(:,:,i) ...
                        + (1-alpha)*Psi_k(:,:,k)*diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))*Psi_k(:,:,k)';
                    
                    Den = alpha*(1+1/gamma)*G_k(:,Nu*(i-1)+u,i)'*(Kz\G_k(:,Nu*(i-1)+u,i));
                    lambda_curr(u,i,k) = 1/real(Den);
                end
            end


        end
 
    [count, sum( abs((lambda_curr(:) - lambda_past(:))./lambda_curr(:)) > th)/(Nu*Nc*Nsc)];
    end
    
    [P_total, P_max, W, P_each] = Precoder_WB_PA1(SP, G, n, lambda_curr, D);
    

    
% %     %% Checks SINR evolution of each iteration
% %     
% %     lambda = lambda_curr;
% % 
% %     SINR = zeros(Nu, Nc, Nsc);
% % 
% % 
% %     for k = 1:Nsc
% %         Psi_k = kron(W_DFT(k,:), eye(Nr));
% % 
% %         G_k = G{k};
% %         lambda_k = lambda(:,:,k);
% %         Lambda_k = diag(lambda_k(:));
% %         for i = 1:Nc
% %             G_i = G_stack(G,i,SP);
% %             % SINR
% %             ICI = G_k(:,:,i)*Lambda_k*G_k(:,:,i)';
% % %             ICI = ICI - lambda_k(u,i)*G_k(:,Nu*(i-1)+u,i)*G_k(:,Nu*(i-1)+u,i)';
% % 
% %             G_incell = G_k(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
% %             lambda_incell = diag(lambda_k(:,i));
% % 
% % 
% % 
% %             F = pinv(alpha^2*ICI + alpha*D(:,:,i) + ... 
% %             alpha*(1-alpha)*Psi_k*(diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))+ kron(eye(Nsc), D(:,:,i)) )*Psi_k')*G_incell; % incell MMSE
% % 
% %             for u = 1:Nu
% %                 SigPow = lambda_incell(u,u)*abs( alpha*F(:,u)'*G_incell(:,u) )^2;
% %                 NoisePow = F(:,u)'* (alpha*D(:,:,i) ...
% %                     + alpha^2*ICI ...
% %                     - alpha^2*G_incell(:,u)*lambda_incell(u,u)*G_incell(:,u)'...
% %                     + alpha*(1-alpha)*Psi_k*(diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))+ kron(eye(Nsc), D(:,:,i)) )*Psi_k') *F(:,u);
% %                 SINR(u,i,k) = SigPow/NoisePow;
% %             end
% %         end
% %     end
% %     
% %     10*log10(mean(abs(SINR(:)))) ;
    %%
    
    D = D_update_WB(SP, D, W);
    [count sum( abs((D(:) - D_past(:))./D(:)) > SP.th2) 10*log10([P_total, P_max])]
end

lambda = lambda_curr;

SINR = zeros(Nu, Nc, Nsc);


for k = 1:Nsc
    Psi = kron(W_DFT, eye(Nr));
    Psi_k = kron(W_DFT(k,:), eye(Nr));
    
    G_k = G{k};
    lambda_k = lambda(:,:,k);
    Lambda_k = diag(lambda_k(:));
    for i = 1:Nc
        G_i = G_stack(G,i,SP);
        % SINR
        ICI = G_k(:,:,i)*Lambda_k*G_k(:,:,i)';
%         ICI = ALL - lambda_k(u,i)*G_k(:,Nu*(i-1)+u,i)*G_k(:,Nu*(i-1)+u,i)';

        G_incell = G_k(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
        lambda_incell = diag(lambda_k(:,i));
        
        
        
        F = pinv(alpha^2*ICI + alpha*D(:,:,i) + ... 
        alpha*(1-alpha)*Psi_k*(diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))+ kron(eye(Nsc), D(:,:,i)) )*Psi_k')*G_incell; % incell MMSE

        for u = 1:Nu
            SigPow = lambda_incell(u,u)*abs( alpha*F(:,u)'*G_incell(:,u) )^2;
            NoisePow = F(:,u)'* (alpha*D(:,:,i) ...
                + alpha^2*ICI ...
                - alpha^2*G_incell(:,u)*lambda_incell(u,u)*G_incell(:,u)'...
                + alpha*(1-alpha)*Psi_k*(diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))+ kron(eye(Nsc), D(:,:,i)) )*Psi_k') *F(:,u);
            SINR(u,i,k) = SigPow/NoisePow;
        end
    end
end


end









