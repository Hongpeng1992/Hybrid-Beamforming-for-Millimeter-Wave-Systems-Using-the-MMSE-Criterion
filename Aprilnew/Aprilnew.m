clear all; close all; clc
% for all
disp(datestr(now));                % show time
%---------structure parameter---------%
SNR_dB = [-24];
SNR = 10 .^ ( SNR_dB/10 );
P = 1;
N_t = 64;                  % Number of transmit antennas
N_r = 64;                  % Number of receive antennas
N_RF = 4;                  % Number of RF chain
Ns = 2;                    % Number of data streams
N_k = 64;                  % the number of subcarrier
fprintf(' simulate the %d * %d * %d %d performance of BER at SNR = %d : %d \n ', N_t,N_r,Ns,N_RF, SNR_dB(1) , SNR_dB(end) );
%%---channel parameter----------%
N_loop =10;               % Number of loop per SNR
cn = 2;
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
        [V_opt, W_opt,MSE] = WBMSEopt( H, Vn ,Ns);     % the optimal
         [V_D(:,:,:,1),V_RF(:,:,1),W_D(:,:,:,1),W_RF(:,:,1),SMSE1(n,:,snr_idx)] = WBEVD(Ns,N_RF,H,Vn,W_opt);
        %  [V_D(:,:,:,1),V_RF(:,:,1),W_D(:,:,:,1),W_RF(:,:,1),SMSE1(n,:,snr_idx)] = WBJUNZhang(N_RF,V_opt,W_opt,Vn,H);
        [V_D(:,:,:,2),V_RF(:,:,2),W_D(:,:,:,2),W_RF(:,:,2),SMSE2(n,:,snr_idx)] = WBMO(Ns,N_RF,H,Vn,W_opt);
        % [V_D(:,:,:,3),V_RF(:,:,3),W_D(:,:,:,3),W_RF(:,:,3),SMSE3(n,:,snr_idx)] =WBOMP(Ns,N_RF,H,Vn,W_opt,AT,AR);
    %    [V_D(:,:,:,3),V_RF(:,:,3),W_D(:,:,:,3),W_RF(:,:,3),SMSE3(n,:,snr_idx)] = NewWBMO(Ns,N_RF,H,Vn,W_opt);
        %% simulate the transmission
        hMod = comm.PSKModulator(4,'BitInput',true,'PhaseOffset',pi/4);
        hDemod = comm.PSKDemodulator('ModulationOrder',4,'BitOutput',true,'PhaseOffset',pi/4);
        data = randi([0 1],N_k*Ns*2,1);
        symbol_t = step(hMod,data);
        symbol_t = reshape(symbol_t,Ns,N_k);
        Noise = sqrt(Vn/2).*(randn(N_r,N_k)+1i*randn(N_r,N_k));
        
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
    BitErr(:,snr_idx) = sum(biterr,2)./(N_k*Ns*2*(N_loop))
end
disp(datestr(now));
