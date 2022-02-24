function [ G ] = Channel_subcarrier(SP, H)

Nsc = SP.Nsc;
tap = SP.tap;

for n = 1:Nsc
    G{n} = zeros(size(H{1}));
    for l = 1:tap
        G{n} = G{n} + H{l}*exp(-1j*2*pi*(n-1)*(l-1)/Nsc);
    end
end
end

