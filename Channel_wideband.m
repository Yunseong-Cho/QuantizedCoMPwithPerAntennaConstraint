function [H] = Channel_wideband( SP )

Nr = SP.Nr; % # Rx Antennas per BS
Nu = SP.Nu;  % # Users per cell
Lp = SP.Lp;
tap = SP.tap;
ch_type = SP.ch_type;
BW = SP.BW;
NF = SP.NF;
dmax = SP.dmax;
dmin = SP.dmin;
Nc = SP.Nc;
Nsc = SP.Nsc;

% SP.Nu = 8; % # Users per cell
% SP.Nr = 8; %ceil(linspace(12,128,15)); % # Antennas to select
% SP.BW = 100*10^6;
% SP.NF = 5; %  Rx Noise Figure in dB
% SP.dmax = 200; % maximum distance in meter
% SP.dmin = 20; % minimum distance in meter
% SP.Lp = 3; % mmWave channel paths
% SP.tap = 4; % chanenl delay taps
% SP.Nsc = 64; % # of subcarrier 
% SP.ch_type = 'mmWave';

switch ch_type
    
    case 'Rayleigh' % sub-6GHz channel
        D = hexcell(SP);
        % large scale fading (pathloss + shadowing)
        for c = 1:Nc
            % large scale fading (pathloss + shadowing)
            d = D(:,:,c); % Distance in meter
            PL = 72.0 + 2.92*10*log10(d) - 8.7*randn(Nu,Nc);
    %         PL = 72.0 + 2.92*10*log10(d) + 8.7*randn(Nms,1);  % Pathloss in dB
            Pnoise = -174 + 10*log10(BW/Nsc) + NF;
            BSant_gain = 15; % sector antenna gain in 3GPP
            gamma = 10.^(-(PL + Pnoise - BSant_gain)/10);
        
        % small scale fading
            for l = 1:tap
                H{l}(:,:,c) = sqrt(1/2)*(randn(Nr, Nu*Nc) + 1j*randn(Nr, Nu*Nc));
                H{l}(:,:,c) = H{l}(:,:,c)*diag(sqrt(gamma(:)));
            end
        end
    case 'mmWave' % mmWave channel
        %D = zeros(Nu,Nc,Nc)
        
        D = hexcell(SP);
        
        for c = 1:Nc
            % large scale fading (pathloss + shadowing)
            d = D(:,:,c); % Distance in meter
            PL = 72.0 + 2.92*10*log10(d) + 8.7*randn(Nu,Nc);
    %         PL = 72.0 + 2.92*10*log10(d) + 8.7*randn(Nms,1);  % Pathloss in dB
            Pnoise = -174 + 10*log10(BW/Nsc) + NF;
            BSant_gain = 15; % sector antenna gain in 3GPP
            gamma = 10.^(-(PL + Pnoise - BSant_gain)/10);

    %         PL = 61.4 + 2*10*log10(d) + 5.7*randn(Nms,1);  % Pathloss in dB
    %         Pnoise = -174 + 10*log10(BW) + NF;
    %         gamma_LOS = 10.^(-(PL + Pnoise)/10);
            
            % small scale fading
            for l = 1:tap
                for c1 = 1:Nc
                    for m = 1:Nu
                        g = 1/sqrt(2)*(randn(Lp,1) + 1j*randn(Lp,1)); % Complex Path Gains
        %                 g(1,1) = g(1,1)*sqrt(gamma_LOS(m));
        %                 g(2:end,1) = g(2:end,1)*sqrt(gamma(m));
                        theta = -1 + 2.*rand(Lp,1); % Angle-of-Arrivals
                        A = SteeringGen(theta, Nr, Lp); % Array Steering Vectors
                        idx = (c1-1)*Nu+m;
                        H{l}(:,idx,c) = sqrt(Nr/Lp)*A*g; % hk
                    end
                end
                H{l}(:,:,c) = H{l}(:,:,c)*diag(sqrt(gamma(:)));
            end
        end
end





