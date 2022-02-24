function [D] = D_update(SP, H, lambda, alpha, D, W)

Nr = SP.Nr;
Nc = SP.Nc;
Nu = SP.Nu;

% LAMBDA = lambda_stack(lambda, SP);


for i = 1:Nc
    
%     grad = zeros(Nr, Nr);
%     
%     for u = 1:Nu
%         grad = grad + W(:,i)*W(:,i)';
%     end
    grad = W(:,:,i)*W(:,:,i)';
    D(:,:,i) = real( D(:,:,i) + SP.mu*diag(diag(grad)) );  
    
    dia = (ones(1,Nr)*diag(D(:,:,i)));
    D(:,:,i) = Nr*D(:,:,i)/dia;

%     while(1) % projection
%         % Step 5: Projection onto hyperplane
% 
%         d = diag(D(:,:,i));
%         one = ones(Nr,1);
%         d = d - (one/norm(one)^2)*(one'*d-Nr);
%         D(:,:,i) = diag(d);
% %         D(:,:,i) = D(:,:,i) - eye(Nr) * max(ones(1,Nr)*diag(D(:,:,i))-Nr,0);
%         %Step 6: Projection onto PSD
%         [V,ev] = eig(D(:,:,i));
%         D(:,:,i) = V*max(ev, 0)*V';
%         if all(diag(ev) >= 0)
%             break;
%         end       
%     end


end


end

