function [ lambda, count , SINR ] = Algo_joint_wb(SP, G, gamma, initlambda, res)
% actual


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


pu = SP.pu;
Nr = SP.Nr;
Nc = SP.Nc;
Nu = SP.Nu;
th = SP.th;
Nsc = SP.Nsc;
algoMax = SP.algoMax;
lambda_curr = initlambda; % initial guess
lambda_past = zeros(Nu,Nc,Nsc);

W_DFT = dftmtx(Nsc)/sqrt(Nsc);

count = 0;
% count2 = 1;
while (sum( abs((lambda_curr(:) - lambda_past(:))./lambda_curr(:)) > th) ~= 0 && count < algoMax)   % relative error, was abs(lambda_c - lambda_p)
%while (1)
    count = count + 1;
    
    lambda_past = lambda_curr;
    
    
    LAMBDA = lambda_stack(lambda_past, SP);
    
    
    k=1;
    i=1;
    u=1;
    
    G_k = G{k}; % channel of k-th subcarrier 
    lambda_k = lambda_past(:,:,k);

    Lambda_k = diag(lambda_k(:));

    Psi = kron(W_DFT, eye(Nr));
    Psi_k = kron(W_DFT(k,:), eye(Nr));

    
    G_i = G_stack(G,i,SP);



    ICI = G_k(:,:,i)*Lambda_k*G_k(:,:,i)';
                      
                
    Kz = alpha*ICI + eye(Nr) + (1-alpha)*Psi_k*diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))*Psi_k';

    D = alpha*(1+1/gamma)*G_k(:,Nu*(i-1)+u,i)'*(Kz\G_k(:,Nu*(i-1)+u,i));
    lambda_curr(u,i,k) = 1/real(D);

  
    
end

%
    lambda = lambda_curr;
    
    Psi = kron(W_DFT, eye(Nr));
    Psi_k = kron(W_DFT(k,:), eye(Nr));
    
    G_k = G{k};
    lambda_k = lambda(:,:,k);
    Lambda_k = diag(lambda_k(:));
    G_i = G_stack(G,i,SP);
        % SINR
    ICI = G_k(:,:,i)*Lambda_k*G_k(:,:,i)';
    G_incell = G_k(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
    lambda_incell = diag(lambda_k(:,i));

    W = pinv(alpha^2*ICI + alpha*eye(Nr) + ... 
    alpha*(1-alpha)*Psi_k*(diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))+eye(Nsc*Nr))*Psi_k')*G_incell; % incell MMSE
    
    SigPow = lambda_incell(u,u)*abs( alpha*W(:,u)'*G_incell(:,u) )^2;
            NoisePow = W(:,u)'* (alpha*eye(Nr) ...
                + alpha^2*ICI ...
                - alpha^2*G_incell(:,u)*lambda_incell(u,u)*G_incell(:,u)'...
                + alpha*(1-alpha)*Psi_k*(diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))+eye(Nsc*Nr))*Psi_k') *W(:,u);
    SINR = SigPow/NoisePow;
    10*log10(abs(SINR));





end









