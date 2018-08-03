clear all; close all; clc
% for all
disp(datestr(now));                % show time
%---------structure parameter---------%
SNR_dB = [-18];
SNR = 10 .^ ( SNR_dB/10 );
P = 1;
N_t = 64;                  % Number of transmit antennas
N_r = 64;                  % Number of receive antennas
N_RF = 2;                  % Number of RF chain
Ns = 2;                    % Number of data streams
N_k = 64;                  % the number of subcarrier
fprintf(' simulate the %d * %d * %d %d performance of BER at SNR = %d : %d \n ', N_t,N_r,Ns,N_RF, SNR_dB(1) , SNR_dB(end) );
%%---channel parameter----------%
N_loop =200;               % Number of loop per SNR
cn = 3;
%%---data matrices----------%
BitErr = zeros(4,length(SNR_dB));     %BER performance
SMSE1 = zeros(N_loop,10,length(SNR_dB));
SMSE2 = zeros(N_loop,10,length(SNR_dB));
SMSE3 = zeros(N_loop,10,length(SNR_dB));
for snr_idx = 1: length(SNR_dB)
    %%---data matrices----------%
    
    biterr = zeros(4 ,N_loop);
    Vn = 1 / 10^(SNR_dB(snr_idx)/10);      % Noise Power
    t1 = clock;
    for n = 1 : N_loop
        %% 先产生频域信道
        [H,AT,AR] = OMPHWB(N_t,N_r);
        %% obtain beamforming matrix
       [V_opt, W_opt,MSE] = WBMSEopt( H, Vn ,Ns);
        [~] = TWBMSEopt( H, Vn ,Ns);     % the optimal
  
       
        [T(1,n)] = TWBEVD(Ns,N_RF,H,Vn,W_opt);
       
       
        %  [V_D(:,:,:,1),V_RF(:,:,1),W_D(:,:,:,1),W_RF(:,:,1),SMSE1(n,:,snr_idx)] = WBJUNZhang(N_RF,V_opt,W_opt,Vn,H);
        t1 = clock;
        %[T(2,n)] = TWBMO(Ns,N_RF,H,Vn,W_opt);
        t2 = clock;
      
        t1 = clock;
        [T(3,n)] =TWBOMP(Ns,N_RF,H,Vn,W_opt,AT,AR);
        t2 = clock;
      
        t1 = clock;
         [~] =YUWEI(Ns,N_RF,H,Vn);
         t2 = clock;
       
        %    [V_D(:,:,:,3),V_RF(:,:,3),W_D(:,:,:,3),W_RF(:,:,3),SMSE3(n,:,snr_idx)] = NewWBMO(Ns,N_RF,H,Vn,W_opt);
        %% simulate the transmission        
      
        %% display
        if (n==10)
            mytoc(N_loop,t1);
        end
    end
    X = sum(T,2)/N_loop
end
disp(datestr(now));
