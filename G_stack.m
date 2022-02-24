function [G_i] = G_stack(G,i,SP)

%i is the present cell index

Nsc = SP.Nsc;
Nc = SP.Nc;
Nu = SP.Nu;

G_i_j = {};
G_i = [];



for j  = 1:Nc
    
    j_idx = 1+(j-1)*Nu:j*Nu;

    for k = 1:Nsc
        g = G{k}; % k-th subcarrier

        G_i_j{k} = g(:,j_idx,i); % cell i -> users in cell j

    end

    % each G_i_j{:}: Nr X Nu
    % blkdiag(G_i_j{:}) : (Nr X Nu) * Nsc
    
    g_i_j = blkdiag(G_i_j{:}); %blkdiag(G_i_j{0}, ..., G_i_j{K-1})

    G_i = [G_i, g_i_j];

end

end

%output = G_i = [g_i_1, ..., g_i_Nc]




