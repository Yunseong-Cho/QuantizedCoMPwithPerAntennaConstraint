function [D] = D_update_WB(SP, D, W)

Nr = SP.Nr;
Nc = SP.Nc;
Nu = SP.Nu;
Nsc = SP.Nsc;



for i = 1:Nc
    
    W_tmp = reshape(W(:,:,i,:), Nr, []);
    grad = diag(diag(W_tmp*W_tmp'));
    
%     for k = 1:Nsc
%         grad = grad + W(:,:,i,k)*W(:,:,i,k)';
%     end

    D(:,:,i) = real( D(:,:,i) + SP.mu*grad );  
    
    dia = (ones(1,Nr)*diag(D(:,:,i)));
    D(:,:,i) = Nr*D(:,:,i)/dia;    
    
%     while(1) % projection
%         % Step 5: Projection onto hyperplane
% 
%         d = diag(D(:,:,i));
%         one = ones(Nr,1);
%         d = d - (one/norm(one)^2)*(one'*d-Nr);
%         D(:,:,i) = diag(d);
%         %Step 6: Projection onto PSD
%         sum(d);
%         D(:,:,i) = max(D(:,:,i), 0);
%         if all(diag(D(:,:,i)) >= 0)
%             break;
%         end 
%     end
    
    
    
end



end

