clear all; close all; clc
% for all
disp(datestr(now));                % show time
%---------structure parameter---------%
SNR_dB = [-10];
P = 1;
N_t = 64;                  % Number of transmit antennas
N_r = 64;                  % Number of receive antennas
N_rf = [7];                  % Number of RF chain
Ns = 4;                    % Number of data streams
N_k = 64;                  % the number of subcarrier
Vn = 1 / 10^(SNR_dB/10); 
%fprintf(' simulate the %d * %d * %d %d performance of BER at SNR = %d : %d \n ', N_t,N_r,Ns,N_RF, SNR_dB(1) , SNR_dB(end) );
%%---channel parameter----------%
N_loop =400;               % Number of loop per SNR
cn = 3;
%%---data matrices----------%
  BitErr = zeros(4,length(N_rf));     %BER performance
  SMSE1 = zeros(N_loop,10,length(N_rf));
    SMSE2 = zeros(N_loop,10,length(N_rf));
    SMSE3 = zeros(N_loop,10,length(N_rf));
for rf_idx = 1: length(N_rf)
    %%---data matrices----------%
    N_RF = N_rf(rf_idx);
    biterr = zeros(4 ,N_loop);
    t1 = clock;
    clear V_D V_RF W_D W_RF
    for n = 1 : N_loop
        %% 先产生频域信道
        [H,AT,AR] = OMPHWB(N_t,N_r);
        %% obtain beamforming matrix
        [V_opt, W_opt] = WBMSEopt( H, Vn ,Ns);     % the optimal
        [V_D(:,:,:,1),V_RF(:,:,1),W_D(:,:,:,1),W_RF(:,:,1),SMSE1(n,:,rf_idx)] = WBEVD(Ns,N_RF,H,Vn,W_opt);
        [V_D(:,:,:,2),V_RF(:,:,2),W_D(:,:,:,2),W_RF(:,:,2),SMSE2(n,:,rf_idx)] = WBMO(Ns,N_RF,H,Vn,W_opt);
        [V_D(:,:,:,3),V_RF(:,:,3),W_D(:,:,:,3),W_RF(:,:,3),SMSE3(n,:,rf_idx)] =WBOMP(Ns,N_RF,H,Vn,W_opt,AT,AR);
        %% simulate the transmission
        hMod = comm.PSKModulator(4,'BitInput',true,'PhaseOffset',pi/4);
        hDemod = comm.PSKDemodulator('ModulationOrder',4,'BitOutput',true,'PhaseOffset',pi/4);
        data = randi([0 1],N_k*Ns*2,1);
        symbol_t = step(hMod,data);
        symbol_t = reshape(symbol_t,Ns,N_k);
        Noise = sqrt(Vn/2).*(randn(N_r,N_k)+1i*randn(N_r,N_k));
        
        symbol_r = zeros(Ns,N_k);
        Hequal = zeros(Ns,Ns,N_k);
        for i = 1 : cn + 1
            symbol_r = zeros(Ns,N_k);
            Hequal = zeros(Ns,Ns,N_k);
            if i ~= (cn + 1)
                for k = 1:N_k
                    Hequal(:,:,k)= W_D(:,:,k,i)'*W_RF(:,:,i)'*H(:,:,k)*V_RF(:,:,i)*V_D(:,:,k,i);
                    symbol_r(:,k) = Hequal(:,:,k)*symbol_t(:,k) + W_D(:,:,k,i)'*W_RF(:,:,i)'*Noise(:,k);
                end
            else
                for k = 1:N_k
                    symbol_r(:,k) = W_opt(:,:,k)'*H(:,:,k)*V_opt(:,:,k)*symbol_t(:,k)+W_opt(:,:,k)'*Noise(:,k);
                end
            end
            symbol_r = step(hDemod,(symbol_r(:)));
            biterr(i,n) = length(find(data~=symbol_r));
        end
        
        %% display
        if (n==10)
            mytoc(N_loop,t1);
        end
    end
    BitErr(:,rf_idx) = sum(biterr,2)./(N_k*Ns*2*(N_loop))
    save RF-1056
end
disp(datestr(now));
