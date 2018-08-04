clear all; close all; clc   
disp(datestr(now));                % show time
%---------structure parameter---------%
SNR_dB = [-24:2:-10];
SNR = 10 .^ ( SNR_dB/10 );
P = 1;
N_t = 64;                  % Number of transmit antennas
N_r = 64;                  % Number of receive antennas
N_RF = 2;                  % Number of RF chain
Ns = 2;                    % Number of data streams
N_s = 64;                  % Number of symbols for one channel realization
fprintf('simulate the %d * %d * %d performance of BER at SNR = %d : %d \n ', N_t,N_r,N_RF, SNR_dB(1) , SNR_dB(end) );
%%---channel parameter----------%
N_loop = 500;          % Number of loop per SNR
cn = 3;                    %comparing varieties
%%---data matrices----------%
BitErr = zeros(cn+1,length(SNR_dB));     %BER performance
CAP = zeros(cn+1,length(SNR_dB));
Jiankong = zeros(1,length(SNR_dB));   % optimal works in partial channel
Diedaicishu = zeros(2,length(SNR_dB)); 
for snr_idx = 1: length(SNR_dB)
    %%---data matrices----------%
    stanMSE = zeros(cn+1 ,N_loop);            %record the MSE   
     Capacity = zeros(cn+1 ,N_loop);
    biterr = zeros(cn+1 ,N_loop);
    dd1 = 0;
    dd2 = 0;
    Vn = 1 / 10^(SNR_dB(snr_idx)/10);      % Noise Power
    t1 = clock;
    n = 1;
    SMSE1 = zeros(N_loop,11);
     SMSE = zeros(N_loop,11);
          SMSE2 = zeros(N_loop,11);
    %MSE1 = zeros(N_loop,20);
    %MSE = zeros(N_loop,20);
    for n = 1 : N_loop        
        %% 先产生频域信道
         H  = Hybrid_channelrlz(N_t,N_r);  
         H = H(:,:,1);
        %% obtain beamforming matrix
        [V_opt, W_opt, MSE_opt] = MSEopt( H, Vn ,Ns);     % the optimal 
         if ( norm(V_opt(:,end))==0 )
             Jiankong(snr_idx) =    Jiankong(snr_idx) + 1;              %avoid the bad case
             continue;
         end;
        [V_D(:,:,1) , V_RF(:,:,1),W_D(:,:,1), W_RF(:,:,1),SMSE1(n,:),n1] =GEini(Ns,N_RF,H,Vn,W_opt);
          %SMSE1 = SMSE1 + MSE1;
%           C1 =  log2(det(eye(Ns) + 1/Ns/Vn * pinv(W_opt) * H * V_opt * V_opt' * H' * W_opt));
%         dd1 = dd1 + n1;  %the number of iteration
         [U,S,V] = svd(H(:,:));
             VC_opt(:,:) = V([1:N_t],[1:Ns]);
             VC_opt(:,:)  = VC_opt(:,:) /norm(VC_opt(:,:) ,'fro');
             WC_opt(:,:) = U([1:N_r],[1:Ns]);
             Capacity(cn+1,n) =  log2(det(eye(Ns) + 1/Ns/Vn * pinv(WC_opt) * H * VC_opt * VC_opt' * H' * WC_opt));
        [V_D(:,:,2) , V_RF(:,:,2),W_D(:,:,2), W_RF(:,:,2)] = JUNZhang(N_RF,V_opt,W_opt); 
         [V_DP,V_RFP,W_DP, W_RFP] = JZPE(N_RF,V_opt,W_opt); 
       %[V_D(:,:,2) , V_RF(:,:,2),W_D(:,:,2), W_RF(:,:,2),SMSE2(n,:),n3] = MOMSEran(Ns,N_RF,H,Vn,W_opt, W_RFP,V_RFP); 
        % [V_D(:,:,3) , V_RF(:,:,3),W_D(:,:,3), W_RF(:,:,3)] = JZPE(N_RF,V_opt,W_opt); 
      
        [V_D(:,:,3) , V_RF(:,:,3),W_D(:,:,3), W_RF(:,:,3),SMSE(n,:),n2] = MOMSEini(Ns,N_RF,H,Vn,W_opt, W_RFP,V_RFP);
        %SMSE = SMSE + MSE;
        %% simulate the transmission
        hMod = comm.PSKModulator(4,'BitInput',true,'PhaseOffset',pi/4);
        hDemod = comm.PSKDemodulator('ModulationOrder',4,'BitOutput',true,'PhaseOffset',pi/4);
        data = randi([0 1],N_s*Ns*2,1);
        symbol_t = step(hMod,data);
        symbol_t = reshape(symbol_t,Ns,N_s);
        Noise = sqrt(Vn/2).*(randn(N_r,N_s)+1i*randn(N_r,N_s));
     
        for i = 1 : cn + 1
            if i ~= cn + 1
                Hequal= W_D(:,:,i)'*W_RF(:,:,i)'*H*V_RF(:,:,i)*V_D(:,:,i);
                symbol_r = Hequal*symbol_t+W_D(:,:,i)'*W_RF(:,:,i)'*Noise;
                stanMSE(i,n) =trace((Hequal - eye(Ns)) * (Hequal - eye(Ns))' + W_D(:,:,i)'*W_RF(:,:,i)'* Vn *eye(N_r)*W_RF(:,:,i)*W_D(:,:,i));
                Capacity(i,n) =  log2(det(eye(Ns) + 1/Ns/Vn * pinv(W_RF(:,:,i)*W_D(:,:,i)) * H * V_RF(:,:,i) * V_D(:,:,i) * V_D(:,:,i)' * V_RF(:,:,i)' * H' *(W_RF(:,:,i)*W_D(:,:,i))));
            else
                Hequal = W_opt'*H*V_opt;
                symbol_r = W_opt'*H*V_opt*symbol_t+W_opt'*Noise;
                stanMSE(i,n) = MSE_opt;
            end
       
        symbol_r = step(hDemod,(symbol_r(:)));
        biterr(i,n) = length(find(data~=symbol_r));
        end
       %% display
        if (mod(n,50) == 0)
           fprintf('已循环 %d 次,error：%d, %d, %d,%d \n',n,sum(biterr(1,:)),sum(biterr(2,:)),sum(biterr(3,:)),sum(biterr(4,:)))
        end
        if (n==10)
          mytoc(N_loop,t1);
        end
    end
    Diedaicishu(1,snr_idx) = dd1/N_loop;
    Diedaicishu(2,snr_idx) = dd2/N_loop;
    BitErr(:,snr_idx) = sum(biterr,2)./(N_s*N_RF*2*(N_loop-Jiankong(snr_idx)))
    CAP(:,snr_idx)=sum(Capacity,2)/N_loop
end
    disp(datestr(now));
