function [P_total, P_max, W, P_each] = Precoder_WB_PA( SP, G, n, lambda, D)


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
Nsc = SP.Nsc;
W_DFT = dftmtx(Nsc)/sqrt(Nsc);
Psi = kron(W_DFT, eye(Nr));


F = zeros(Nr, Nu, Nc, Nsc);
W = zeros(Nr, Nu, Nc, Nsc);

LAMBDA = lambda_stack(lambda, SP);


%% To get combiner in UL

for k = 1:Nsc
    G_k = G{k};
    lambda_k = lambda(:,:,k);
    Lambda_k = diag(lambda_k(:));
    Psi_k = kron(W_DFT(k,:), eye(Nr));
    
    for i = 1:Nc
        G_i = G_stack(G,i,SP);

        for u = 1:Nu

    %         ICI = H(:,:,i)*diag(lambda(:))*H(:,:,i)';
            ALL = G_k(:,:,i)*Lambda_k*G_k(:,:,i)';
            ICI = ALL - lambda_k(u,i)*G_k(:,Nu*(i-1)+u,i)*G_k(:,Nu*(i-1)+u,i)';
            G_incell = G_k(:,Nu*(i-1)+u,i);
%             G_incell = G_k(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
            F(:, u , i, k) = pinv(alpha^2*ICI ...
                                + alpha^2*D(:,:,i) ...
                                + alpha*(1-alpha)*Psi_k*(diag(diag(Psi'*G_i*LAMBDA*G_i'*Psi))+ kron(eye(Nsc), D(:,:,i)))*Psi_k')*G_incell; % incell MMSE
        end

    end
end


%% Corollary 3, DL  
        Sigma = [];
        
        for k = 1:Nsc
            Sigma_k = zeros(Nu*Nc, Nu*Nc);
            G_k = G{k};
            for i = 1:Nc
                for j = 1:Nc
                    %Each Sigma matrix

                    for u = 1:Nu
                        for v = 1:Nu

                            row = (i-1)*Nu+u;
                            col = (j-1)*Nu+v;


                            if i==j && u==v
                                f = F(:,u,i,k);
                                g = G_k(:,Nu*(i-1)+u,i);
                                Sigma_k(row,col) = (alpha^2/gamma)*(abs(f'*g))^2;
                                for n = 1:Nsc
                                    fn = F(:,u,i,n);
                                    Psi_n = kron(W_DFT(n,:), eye(Nr));
                                    gg = [zeros((n-1)*Nr,1); g; zeros((Nsc-n)*Nr,1)];
                                    Sigma_k(row,col) = Sigma_k(row,col) - alpha*(1-alpha)*fn'*Psi_n*diag(diag(Psi'*gg*gg'*Psi))*Psi_n'*fn;      
                                end
                                
                            else
                                f = F(:,v,j,k);
                                g = G_k(:,Nu*(i-1)+u,j);
                                Sigma_k(row,col) = -(alpha^2)*(abs(f'*g))^2;  
                                for n = 1:Nsc
                                    fn = F(:,v,j,n);
                                    Psi_n = kron(W_DFT(n,:), eye(Nr));
                                    gg = [zeros((n-1)*Nr,1); g; zeros((Nsc-n)*Nr,1)];
                                    Sigma_k(row,col) = Sigma_k(row,col) - alpha*(1-alpha)*fn'*Psi_n*diag(diag(Psi'*gg*gg'*Psi))*Psi_n'*fn;      
                                end
                            end
                        end
                    end
                end
            end
            if k == 1
                Sigma = Sigma_k;
            else
                Sigma = [Sigma zeros((k-1)*Nu*Nc,Nu*Nc);zeros(Nu*Nc,(k-1)*Nu*Nc) Sigma_k];
            end
                
        end
        
        tau = pinv(Sigma)*ones(Nsc*Nu*Nc,1);
        
        for k = 1:Nsc
            for i = 1:Nc
                for u = 1:Nu
                    idx = (k-1)*Nc+(i-1)*Nu+u;
                    W(:,u,i,k) = sqrt(tau(idx))*F(:,u,i,k);
                end
            end
        end
%  
    
    
    P_total = W(:)' * W(:);
    
    P_each = [];
    
    for i = 1:Nc
        P=0;
        for k = 1:Nsc
            P = P + (alpha/Nsc)*W(:,:,i,k)*W(:,:,i,k)';
        end
        
        P_each = [P_each; diag(P)];
    end
    
    P_max = max(P_each);

    
    
%     SINR = zeros(Nu,Nc);
%     for i = 1:Nc
% 
%         % SINR
%         ICI = H(:,:,i)*diag(lambda(:))*H(:,:,i)';
%         H_incell = H(:,Nu*(i-1)+1:Nu*(i-1)+Nu,i);
%         lambda_incell = diag(lambda(:,i));
%         W = pinv(alpha^2*ICI + alpha*eye(Nr) + alpha*(1-alpha)*diag(diag(ICI)))*H_incell; % incell MMSE
% 
%         for u = 1:Nu
%             SigPow = lambda_incell(u,u)*abs( alpha*W(:,u)'*H_incell(:,u) )^2;
%             NoisePow = W(:,u)'* (alpha*eye(Nr) ...
%                 + alpha^2*ICI - alpha^2*H_incell(:,u)*lambda_incell(u,u)*H_incell(:,u)'...
%                 + alpha*(1-alpha)*diag(diag(H(:,:,i)*diag(lambda(:))*H(:,:,i)'))) *W(:,u);
%             SINR(u,i) = SigPow/NoisePow;
%         end
%     end

    

end














