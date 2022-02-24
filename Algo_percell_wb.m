function [ lambda, cnt, SINR ] = Algo_percell_wb(SP, G, gamma, initlambda, res)

switch res
    case 'low'
        b = SP.b;
        bTable = [0.3634, 0.1175, 0.03454, 0.009497, 0.002499];
        if b > 5
            alpha = (1 - pi*sqrt(3)/2*2.^(-2*b));
        else
            alpha = (1 - bTable(floor(b)));
        end
        
    case 'infinite'
        alpha = 1;
end

Nr = SP.Nr;
Nc = SP.Nc;
Nu = SP.Nu;
th = SP.th;
Nsc = SP.Nsc;
algoMax = SP.algoMax;
algoMax2 = SP.algoMax2;

lambda_curr = initlambda; % initial guess
lambda_past = zeros(Nu,Nc,Nsc);

cnt = 0;
count = 0;
count2 = 0;

W_DFT = dftmtx(Nsc)/sqrt(Nsc);


while (sum( abs((lambda_curr(:) - lambda_past(:))./lambda_curr(:)) > th) ~= 0 && count2 < algoMax ) % all cells  % relative error, was abs(lambda_c - lambda_p)
    
    lambda_past = lambda_curr;
    LAMBDA = lambda_stack(lambda_past, SP); % diag over all powers
    
    for i = 1:Nc
        for k = 1:Nsc
            
            % channel
            G_k = G{k}; % channel of k-th subcarrier, Nb X Nu*Nc X Nc
            lambda_k = lambda_past(:,:,k);

            Psi = kron(W_DFT, eye(Nr));
            Psi_k = kron(W_DFT(k,:), eye(Nr));    
            G_i = G_stack(G,i,SP); % Nsc*Nr X Nsc*Nu*(Nc) 
            
            % Noise estimate
            lambda_out = lambda_past(:,:,k);
            lambda_out = lambda_out(:);
            lambda_out(Nu*(i-1)+1:Nu*(i-1)+Nu) = 0;
 
            noise = diag(alpha^2*G_k(:,:,i)*diag(lambda_out)*G_k(:,:,i)') + alpha*eye(Nr); %replace G_i with 16*12
            noise_q = alpha*(1-alpha)*Psi_k*diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))*Psi_k';

            lambda_k_cell_c = lambda_k(:,i); % Nu X 1
            lambda_k_cell_p = zeros(Nu,1);

            while (sum( abs((lambda_k_cell_c - lambda_k_cell_p)./lambda_k_cell_c) > th) ~= 0)  
                lambda_k_cell_p = lambda_k_cell_c;
                
                lambda = lambda_past(:,:,k); % the values generated before getting into the inner loop
                lambda(:,i) = lambda_k_cell_p; % Replace the value of corresponding cell only
                Lambda_percell = diag(lambda(:));

                ICI = G_k(:,:,i)*Lambda_percell*G_k(:,:,i)';
                G_incell = G_k(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
                lambda_incell = diag(lambda_k_cell_p);
                
                W = pinv(alpha^2*G_incell*lambda_incell*G_incell' + diag(noise + noise_q))*G_incell; % incell MMSE

                for u = 1:Nu
                    SigPow = abs( alpha*W(:,u)'*G_incell(:,u) )^2;
                    NoisePow = W(:,u)'* (alpha*eye(Nr) ...
                        + alpha^2*ICI ...
                        - alpha^2*G_incell(:,u)*lambda_incell(u,u)*G_incell(:,u)'...
                        + alpha*(1-alpha)*diag(diag(G_k(:,:,i)*Lambda_percell*G_k(:,:,i)')))*W(:,u);
                    lambda_k_cell_c(u) = gamma*real(NoisePow/SigPow);

                end
                

            end
            
        lambda_curr(:,i,k) = lambda_k_cell_c;
    
        end
        count = count + 1;
    end
    count2 = count2 + 1;
end

%% Compute SINR


lambda = lambda_curr;

SINR = zeros(Nu,Nc, Nsc);
LAMBDA = lambda_stack(lambda, SP); % diag over all powers

for k = 1:Nsc
    G_k = G{k}; % channel of k-th subcarrier, Nb X Nu*Nc X Nc
    Psi = kron(W_DFT, eye(Nr));
    Psi_k = kron(W_DFT(k,:), eye(Nr));    
    G_i = G_stack(G,i,SP); % Nsc*Nr X Nsc*Nu*(Nc) 
    
    lambda_k = lambda(:,:,k);
    Lambda_k = diag(lambda_k(:));

    for i = 1:Nc
        % Noise estimate
        lambda_out = lambda(:,:,k);
        lambda_out(Nu*(i-1)+1:Nu*(i-1)+Nu) = 0;
        noise = diag(alpha^2*G_k(:,:,i)*diag(lambda_out(:))*G_k(:,:,i)') + alpha*eye(Nr); %replace G_i with 16*12
        noise_q = alpha*(1-alpha)*Psi_k*diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))*Psi_k';

        % SINR
        G_incell = G_k(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
        lambda_incell = diag(lambda_k(:,i));
        
        W = pinv(alpha^2*G_incell*lambda_incell*G_incell' + diag(noise + noise_q))*G_incell; % incell MMSE
        ICI = G_k(:,:,i)*Lambda_k*G_k(:,:,i)';
        for u = 1:Nu
            SigPow = lambda_incell(u,u)*abs( alpha*W(:,u)'*G_incell(:,u) )^2;
            NoisePow = W(:,u)'* (alpha*eye(Nr) ...
                + alpha^2*ICI ...
                - alpha^2*G_incell(:,u)*lambda_incell(u,u)*G_incell(:,u)'...
                + alpha*(1-alpha)*diag(diag(G_k(:,:,i)*Lambda_percell*G_k(:,:,i)'))) *W(:,u);
            SINR(u,i,k) = SigPow/NoisePow;
        end
    end
end
%10*log10(mean(abs(SINR(:))))

end



