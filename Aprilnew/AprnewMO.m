%clear all; close all; clc
% for WBMO
disp(datestr(now));                % show time
%---------structure parameter---------%
SNR_dB = [-10];
SNR = 10 .^ ( SNR_dB/10 );
P = 1;
N_t = 64;                  % Number of transmit antennas
N_r = 64;                  % Number of receive antennas
N_RF = 2;                  % Number of RF chain
Ns = 2;                    % Number of data streams
N_k = 64;                  % the number of subcarrier
fprintf('WBMO simulate the %d * %d * %d %d performance of BER at SNR = %d : %d \n ', N_t,N_r,Ns,N_RF, SNR_dB(1) , SNR_dB(end) );
%%---channel parameter----------%
N_loop =1000;               % Number of loop per SNR
%%---data matrices----------%
BitErr = zeros(1,length(SNR_dB));     %BER performance
%load('HZmat.mat')
for snr_idx = 1: length(SNR_dB)
    %%---data matrices----------%
    stanMSE = zeros(1 ,N_loop);            %record the MSE
    SMSE1 = zeros(N_loop,10,length(SNR_dB));
    biterr = zeros(1 ,N_loop);
    Vn = 1 / 10^(SNR_dB(snr_idx)/10);      % Noise Power
    t1 = clock;
    for n = 1 : N_loop
        %% 先产生频域信道
        x = mod((n-1),1000)+1;
        H = HZ(:,:,:,x);
        %% obtain beamforming matrix
        [V_opt, W_opt] = WBMSEopt( H, Vn ,Ns);     % the optimal
        [V_D,V_RF,W_D,W_RF,SMSE1(n,:,snr_idx)] = WBMO(Ns,N_RF,H,Vn,W_opt);
        %% simulate the transmission
        hMod = comm.PSKModulator(4,'BitInput',true,'PhaseOffset',pi/4);
        hDemod = comm.PSKDemodulator('ModulationOrder',4,'BitOutput',true,'PhaseOffset',pi/4);
        data = randi([0 1],N_k*Ns*2,1);
        symbol_t = step(hMod,data);
        symbol_t = reshape(symbol_t,Ns,N_k);
        Noise = sqrt(Vn/2).*(randn(N_r,N_k)+1i*randn(N_r,N_k));
        
        symbol_r = zeros(Ns,N_k);
        Hequal = zeros(Ns,Ns,N_k);
        for k = 1:N_k
            Hequal(:,:,k)= W_D(:,:,k)'*W_RF'*H(:,:,k)*V_RF*V_D(:,:,k);
            symbol_r(:,k) = Hequal(:,:,k)*symbol_t(:,k) + W_D(:,:,k)'*W_RF'*Noise(:,k);
        end
        
        symbol_r = step(hDemod,(symbol_r(:)));
        biterr(n) = length(find(data~=symbol_r));
        
        %% display
        if (n==10)
            mytoc(N_loop,t1);
        end
    end
    BitErr(snr_idx) = sum(biterr)./(N_k*Ns*2*(N_loop))
end
disp(datestr(now));
