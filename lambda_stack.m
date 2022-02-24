function LAMBDA = lambda_stack(lambda, SP);


Nsc = SP.Nsc;
Nc = SP.Nc;

lambda_i = {};
Lambda_i = {};

for i = 1:Nc

    for k = 1:Nsc  
        
        lambda_i{k} = diag(lambda(:,i,k));
        
    end
    
    Lambda_i{i} = blkdiag(lambda_i{:});
    
    
    
end


LAMBDA = blkdiag(Lambda_i{:});



