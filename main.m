clear all;
close all
% joint only
rng(3)
SP.Nr = 32;  Nr = SP.Nr;     % # Rx Antennas per BS
SP.Nc = 5;  Nc = SP.Nc;      % # Cells
SP.Nu = 2;  Nu = SP.Nu;      % # Users per cell
SP.Nsc = 1; Nsc = SP.Nsc;

% 
%42.5402      16.3631 vs 41.932       13.506  (noPA vs PA at -5 db) inf...
%42.833       16.503  vs 42.228       13.568  (noPA vs PA at -5 db) 3 bits...

%39.85       13.543  vs 39.255       10.601 (at -8 dB) 3bits
 

% % % PA(:) = SINR_PA{:}
% % % PA(:) = 10.^((10*log10(PA(:))+0.16)/10)
% % % cdfplot(10*log10(PA(:)))
% 
% hold on
% PA(:) = SINR_PA{:}
% cdfplot(10*log10(PA(:)))


SP.b = 3;       % # Quantization Bits
SP.pu_dBm = 0;
SP.pu = 10.^(SP.pu_dBm/10); % mW 


% SP.Nbs = 32; %ceil(linspace(16,128,10)); % # Rx Antennas per BS
% SP.Nms = 8; % # Users per cell
% SP.Nr = 8; %ceil(linspace(12,128,15)); % # Antennas to select
% SP.b = 3; % # Quantization Bits
% SP.p_dBm = linspace(10,40,15?mW
% SP.BW = 100*10^6;
% SP.NF = 5; %  Rx Noise Figure in dB
% SP.dmax = 200; % maximum distance in meter
% SP.dmin = 20; % minimum distance in meter
% SP.Lp = 3; % mmWave channel paths
% SP.tap = 4; % chanenl delay taps
% SP.Nsc = 64; % # of subcarrier 
% SP.Nmcmc = 6; % # of MCMC samples 
% SP.tmax = 3; % # of trials
% SP.tau = 1; % Rate constant (need to use simulated annealing)
% SP.ch_type = 'mmWave';
% SP.iterMax = 100;
% 

dB_min = 12;
dB_max = 12;
SP.gamma_dB = linspace(dB_min, dB_max, 1);
SP.gamma_dB = [2];
SP.gamma = 10.^(SP.gamma_dB/10);

% SP.gamma = linspace(10^(dB_min/10),  10^(dB_max/10), 4);
% 
% SP.gamma_dB = 10*log10(SP.gamma);


SP.th = 1e-2;              % Iteration Threshold
SP.th2 = 1e-3;              % Threshold for D update
SP.TcMax = 1;
SP.initpower  = 10;
SP.algoMax = 500;
SP.algoMax2 = 2500;
SP.iterMax = 1;

% Channel Parameters
SP.BW = 10*10^6;
SP.NF = 5;
SP.dmax = 200;
SP.dmin = 50;
SP.d0 = 100;
SP.rho = 3;
SP.s = 8.7;
SP.LD = 4*10^8/(2.4*10^9);
SP.a = 0.8; % First-order Gauss-Markov model for small scale fading
SP.Lp = 3;  
SP.tap = 3;
SP.ch_type = 'Rayleigh';

SP.mu = 5e-1; % step size
%-0.72 w/ 5e-6
%-0.71 w/ 1e-6

%-0.6693 w/ 4 4 1e-6

%%

ADC = [3];


P_joint_perantenna = zeros(length(SP.gamma),length(ADC));
SINR_joint_perantenna = zeros(length(SP.gamma),length(ADC));

P_joint_no_perantenna = zeros(length(SP.gamma),length(ADC));
SINR_joint_no_perantenna = zeros(length(SP.gamma),length(ADC));

P_percell_result = zeros(length(SP.gamma),length(ADC));
SINR_percell_result = zeros(length(SP.gamma),length(ADC));

P_max_PA = zeros(length(SP.gamma),length(ADC));
P_max_noPA = zeros(length(SP.gamma),length(ADC));

PAPR_noPA = zeros(length(SP.gamma),length(ADC));
PAPR_PA = zeros(length(SP.gamma),length(ADC));
%%

f = waitbar(0);
now = 0;
all = length(SP.gamma)*SP.iterMax*(length(ADC)+1);

%rng(9)
P_noPA = {};
P_PA = {};
SINR_noPA = {};
SINR_PA = {};

for iter = 1:SP.iterMax 
    iter
    [H] = Channel_wideband( SP );
    [G] = Channel_subcarrier(SP, H);
    
    for n = 1:length(SP.gamma)   
        SP.initpower = 10^((29+SP.gamma_dB(n))/10)/SP.Nu/SP.Nc/SP.Nsc;
        initlambda = SP.initpower*ones(SP.Nu,SP.Nc,SP.Nsc);

       

        now = now+1;
        waitbar(now/all,f,num2str(now*100/all));
        
        for bb = 1:length(ADC)
            
            SP.b = ADC(bb);
            
            %%
            SP.perantenna = 0;
            [P_jl, SINR_jl, ~, lambda] = CoMP_WB(SP, G, n, 'joint',initlambda);
            P_joint_no_perantenna(n, bb) = P_joint_no_perantenna(n, bb) + (P_jl - P_joint_no_perantenna(n, bb))/iter;
            SINR_joint_no_perantenna(n, bb) = SINR_joint_no_perantenna(n, bb)+ (mean(SINR_jl(:)) - SINR_joint_no_perantenna(n, bb))/iter;
            
            [P_DL2, P_max2, ~, P_each2] = Precoder_WB(SP, G, n, lambda);
            
            PAPR_noPA(n, bb) = PAPR_noPA(n, bb) + (max(P_each2(1:Nc))/mean(P_each2(1:Nc)) - PAPR_noPA(n, bb))/iter;
            
            P_noPA{n,iter} = P_each2(:);
            SINR_noPA{n,iter} = SINR_jl(:);

            
            noPA = 10*log10(abs([P_DL2, P_max2]))
            
%             noPA = [0,0];
            
            %%
            SP.perantenna = 1;
            [P_jl_p, SINR_jlp, ~, lambda, D] = CoMP_WB(SP, G, n, 'joint',initlambda, noPA);
            P_joint_perantenna(n, bb) = P_joint_perantenna(n, bb) + (P_jl_p - P_joint_perantenna(n, bb))/iter;
            SINR_joint_perantenna(n, bb) = SINR_joint_perantenna(n, bb)+ (mean(SINR_jlp(:)) - SINR_joint_perantenna(n, bb))/iter;
            
            [P_DL1, P_max1, ~, P_each1] = Precoder_WB_PA1(SP, G, n, lambda, D);
            P_PA{n,iter} = P_each1(:);
            SINR_PA{n,iter} = SINR_jlp(:);
            
            PAPR_PA(n, bb) = PAPR_PA(n, bb) + (max(P_each1(1:Nc))/mean(P_each1(1:Nc)) - PAPR_PA(n, bb))/iter;
            
            
            PA = 10*log10(abs([P_DL1, P_max1]))
                

            
            %%
            P_max_PA(n, bb) = P_max_PA(n, bb) + (P_max1 - P_max_PA(n, bb))/iter;
            P_max_noPA(n, bb) = P_max_noPA(n, bb) + (P_max2 - P_max_noPA(n, bb))/iter;
            
    %         [P_pl, SINR_pl, ~, ~] = CoMP_wb(SP, G, n, 'percell', 'low', initlambda);
    %         P_pl_result(n, bb) = P_pl_result(n, bb) + (P_pl - P_pl_result(n, bb))/iter;
    %         SINR_pl_result(n, bb) = SINR_pl_result(n, bb)+ (mean(SINR_pl(:)) - SINR_pl_result(n, bb))/iter;


            now = now+1;
            waitbar(now/all,f,num2str(now*100/all));            
%             total_transmit_UL = [10*log10(P_jl_p) 10*log10(P_jl)]

            
%             achieved_SINR = [10*log10(SINR_joint_perantenna) 10*log10(SINR_joint_no_perantenna)]

       
        
        end
        


        
    end    
end

10*log10([PAPR_PA, PAPR_noPA])

%%
figure;
plot(SP.gamma_dB, 10*log10(P_max_noPA))

hold on
plot(SP.gamma_dB, 10*log10(P_max_PA))

legend('noPA','PA')
xlabel('target SINR [dB]')
ylabel('Max antenna power')


figure;
for n = 1:length(SP.gamma)  
    subplot(1, length(SP.gamma), n)
    noPA = SINR_noPA{n,:};
    PA = SINR_PA{n,:};
    
    cdfplot(10*log10(noPA(:)))
    hold on
    cdfplot(10*log10(PA(:)))
    legend('w/o PA','w/ PA')
    axis([SP.gamma_dB(n)-5, SP.gamma_dB(n)+5, -inf, inf])
end

figure;

plot(SP.gamma_dB, PAPR_noPA)
hold on
plot(SP.gamma_dB, PAPR_PA)
legend('noPA','PA')
xlabel('target SINR [dB]')
ylabel('PAPR')


figure;
cdfplot(10*log10(abs(P_each2))) % no PA
hold on
cdfplot(10*log10(abs(P_each1))) % PA
legend('no PA','PA')
% P_p = 10*log10([P_pl_result, P_pi_result])
% P_j = 10*log10([P_jl_result, P_ji_result])

%%
% figure
% hold on
% % plot(SP.gamma_dB, P_pl_result(:,1), 'b:')
% % plot(SP.gamma_dB, P_pl_result(:,2), 'b--')
% % plot(SP.gamma_dB, P_pi_result, 'b-')
% plot(SP.gamma_dB, 10*log10(P_joint_perantenna(:,1)), 'bo:')
% plot(SP.gamma_dB, 10*log10(P_joint_perantenna(:,2)), 'bo-.')
% plot(SP.gamma_dB, 10*log10(P_joint_perantenna(:,3)), 'bo--')
% plot(SP.gamma_dB, 10*log10(P_joint_perantenna(:,4)), 'bo-')
% 
% plot(SP.gamma_dB, 10*log10(P_joint_no_perantenna(:,1)), 'ro:')
% plot(SP.gamma_dB, 10*log10(P_joint_no_perantenna(:,2)), 'ro-.')
% plot(SP.gamma_dB, 10*log10(P_joint_no_perantenna(:,3)), 'ro--')
% plot(SP.gamma_dB, 10*log10(P_joint_no_perantenna(:,4)), 'ro-')
% 
% xlabel('Target SINR [dB]')
% ylabel('Total Transmit Power [dBm]')
% title(['(Nr, Nc, Nu, Nsc, tap)=', num2str([SP.Nr, SP.Nc, SP.Nu, SP.Nsc, SP.tap])])
% 
% legend('Q-iCoMP ($b=2$)','Q-iCoMP ($b=3$) ','Q-iCoMP ($b=4$)','Q-iCoMP ($b=\infty$)','Interpreter','latex')
% grid on


