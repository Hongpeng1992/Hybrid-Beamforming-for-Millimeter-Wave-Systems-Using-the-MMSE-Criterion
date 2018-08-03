clear all; close all; clc
% for all
disp(datestr(now));                % show time
%---------structure parameter---------%
SNR_dB = [-15];
SNR = 10 .^ ( SNR_dB/10 );
P = 1;
N_t = 64;                  % Number of transmit antennas
N_r = 64;                  % Number of receive antennas
N_RF = 5;                  % Number of RF chain
Ns = 2;                    % Number of data streams
N_k = 64;                  % the number of subcarrier
fprintf(' simulate the %d * %d * %d %d performance of BER at SNR = %d : %d \n ', N_t,N_r,Ns,N_RF, SNR_dB(1) , SNR_dB(end) );
%%---channel parameter----------%
N_loop =100;               % Number of loop per SNR
cn = 3;
%%---data matrices----------%
BitErr = zeros(cn+1,length(SNR_dB));     %BER performance
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
        load('H.mat')
        %% obtain beamforming matrix
        [V_opt, W_opt] = WBMSEopt( H, Vn ,Ns);     % the optimal
        for i = 1:N_k
            [U,S,V] = svd(H(:,:,i));
            VC_opt(:,:,i) = V([1:N_t],[1:Ns]);
            VC_opt(:,:,i)  = VC_opt(:,:,i) /norm(VC_opt(:,:,i) ,'fro');
            
            WC_opt(:,:,i) = U([1:N_r],[1:Ns]);
            C1(i) =  log2(det(eye(Ns) + 1/Vn * pinv(WC_opt(:,:,i)) * H(:,:,i) * VC_opt(:,:,i) * VC_opt(:,:,i)' * H(:,:,i)' * WC_opt(:,:,i)));
        end
        C(1,n) = sum(C1)/N_k;
        
        [V_D(:,:,:,1),V_RF(:,:,1),W_D(:,:,:,1),W_RF(:,:,1)]= WBOMP(Ns,N_RF,H,Vn,W_opt,AT,AR);
        [V_D(:,:,:,2),V_RF(:,:,2),W_D(:,:,:,2),W_RF(:,:,2)]= WBEVD(Ns,N_RF,H,Vn,W_opt);
        [V_D(:,:,:,3),V_RF(:,:,3),W_D(:,:,:,3),W_RF(:,:,3)] =YUWEI(Ns,N_RF,H,Vn);
       % [V_D(:,:,:,4),V_RF(:,:,4),W_D(:,:,:,4),W_RF(:,:,4)]= WBJUNZhang(N_RF,V_opt,W_opt,Vn,H);
      %  [V_D(:,:,:,5),V_RF(:,:,5),W_D(:,:,:,5),W_RF(:,:,5)] = WBMO(Ns,N_RF,H,Vn,W_opt);
        for m = 1:cn
            for i = 1:N_k
                C2(i) = log2(det(eye(Ns) + 1/Vn * pinv(W_RF(:,:,m)*W_D(:,:,i,m)) * H(:,:,i) * V_RF(:,:,m)*V_D(:,:,i,m) * (V_RF(:,:,m)*V_D(:,:,i,m) )' * H(:,:,i)' *W_RF(:,:,m)*W_D(:,:,i,m)));
                %  C4(i) = log2(det(eye(Ns) + 1/Ns/Vn * pinv(WC_opt(:,:,i)) * H(:,:,i) * V_RF(:,:,3)*V_D(:,:,i,3) * (V_RF(:,:,3)*V_D(:,:,i,3) )' * H(:,:,i)' *WC_opt(:,:,i)));
            end
            C(m+1,n) = sum(C2)/N_k;
        end
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
    % BitErr(:,snr_idx) = sum(biterr,2)./(N_k*Ns*2*(N_loop))
    CCAP(:,snr_idx) = sum(C,2)/N_loop
    BitErr(:,snr_idx) = sum(biterr,2)./(N_k*Ns*2*(N_loop))
end
disp(datestr(now));
