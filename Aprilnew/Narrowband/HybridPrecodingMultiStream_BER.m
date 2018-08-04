clear all;
close all;clc
% rng('default')
%------description------------
%Hybrid MMSE Precoding and combining BER
%Hybrid Precoding at BS
%MultiAnt and single RF chain at MS
%-------END-------------------

% SNR_dB = [ -30 ];
SNR_dB = [ -24];
disp(datestr(now));                % show time
fprintf('simulate the performance of BER at SNR = %d : %d \n ', SNR_dB(1) , SNR_dB(end) ) ;
SNR = 10 .^ ( SNR_dB/10 );
P = 1;                     % total transmit power
%---------structure parameter---------%
N_t = 64;                  % Number of transmit antennas
N_r = 64;                   % Number of receive antennas
N_RFt = 2;                 % Number of transmit RF chain
N_RFr = 2;                 % Number of receive RF chain
N_u = 1;                   % Number of users
N_s = 2;                   % Number of data streams
N_sym = 64;                % Number of symbols
%%---channel parameter----------%
L = 20;
N_loop = 20;            % Number of loop per SNR

hMod = comm.PSKModulator(4,'BitInput',true,'PhaseOffset',pi/4);
hDemod = comm.PSKDemodulator('ModulationOrder',4,'BitOutput',true,'PhaseOffset',pi/4);
% constellation(hDemod)
BitErr = zeros(1,length(SNR_dB));
BitErr_FD = zeros(1,length(SNR_dB));
BitErrY = zeros(1,length(SNR_dB));
BitErrPE = zeros(1,length(SNR_dB));
sumMSE = zeros(1,length(SNR_dB));
sumMSE_FD = zeros(1,length(SNR_dB));
sumMSEY = zeros(1,length(SNR_dB));
sumMSEPE = zeros(1,length(SNR_dB));

for snr_idx = 1: length(SNR_dB)
    biterr = zeros(1,N_loop);
    biterr_FD = zeros(1,N_loop);
    biterrY = zeros(1,N_loop);
    biterrPE = zeros(1,N_loop);
    mse = zeros(1,N_loop);
    mse_FD = zeros(1,N_loop);
    mseY = zeros(1,N_loop);
    msePE = zeros(1,N_loop);
   
    Vn = 1 / 10^(SNR_dB(snr_idx)/10);               % Noise Power
    for n = 1: N_loop
        % 产生信道
       % [ H,W_codebook,F_codebook ] = ChannelULA( N_r,N_t,N_u,L );
        [ H ] = Hybrid_channelrlz(64,64);
         [V_opt, W_opt, MSE_opt] = MSEopt( H, Vn ,N_s); 
        % Fully Digital
       % [ W_FD,V_FD,beta_FD,mse_FD(n) ] = GED_OptReceiver( N_r,N_t,N_u,N_s,Vn,P,H );
        % 提出的方案
%         [ W_RF,V_RF,W_D,V_D,beta,mse(n) ] = GED_Opt( N_r,N_t,N_RFt,N_RFr,N_u,N_s,Vn,P,H );
%         [ W_RF1,V_RF1,W_D1,V_D1,beta1,mseY(n) ] = GED_Opt1( N_r,N_t,N_RFt,N_RFr,N_u,N_s,Vn,P,H );
        [ W_RF,V_RF,W_D,V_D,beta,mse(n) ] = Uni_GED_Opt( N_r,N_t,N_RFt,N_RFr,N_u,N_s,Vn,P,H,W_opt );
        % 对比曲线
        [ W_RF_PE,V_RF_PE,W_D_PE,V_D_PE,beta_PE,msePE(n) ] = PE_GED_Opt( W_FD,V_FD./beta_FD,N_r,N_t,N_RFt,N_RFr,N_u,N_s,Vn,P,H );

        
        %% data transmission processing
        data = randi([0 1],N_s*N_sym*N_u*2,1);
        symbol_t = step(hMod,data);
        symbol_t = reshape(symbol_t,N_s*N_u,N_sym);
        Noise = sqrt(Vn/2).*(randn(N_r*N_u,N_sym)+1i*randn(N_r*N_u,N_sym));
        Heq = ChannelEqu( H );      %用户总的信道
%         % 提出的算法
        symbol_r = Heq*V_RF*V_D*symbol_t+Noise;
        W_RFeq = blkdiagExtend(W_RF);
        W_Deq = blkdiagExtend(W_D);
        desymbol = beta^(-1)*W_Deq'*W_RFeq'*symbol_r;
        desymbol = desymbol(:);
        symbol_r = step(hDemod,desymbol);
        biterr(n) = length(find(data~=symbol_r));
        % 接收端Fully Digital
        symbol_r_FD = Heq*V_FD*symbol_t+Noise;
        W_FD = blkdiagExtend(W_FD);
        desymbol_FD = beta_FD^(-1)*W_FD'*symbol_r_FD;
        desymbol_FD = desymbol_FD(:);
        symbol_r_FD = step(hDemod,desymbol_FD);
        biterr_FD(n) = length(find(data~=symbol_r_FD));
%         % 发送端hybrid接收端单位阵
%         symbol_rY = Heq*V_RF1*V_D1*symbol_t+Noise;
%         W_RFeq1 = blkdiagExtend(W_RF1);
%         W_Deq1 = blkdiagExtend(W_D1);
%         desymbolY = beta1^(-1)*W_Deq1'*W_RFeq1'*symbol_rY;
%         desymbolY = desymbolY(:);
%         symbol_rY = step(hDemod,desymbolY);
%         biterrY(n) = length(find(data~=symbol_rY));
        % 发送端hybrid接收端PE
        symbol_rPE = Heq*V_RF_PE*V_D_PE*symbol_t+Noise;
        W_RFeq_PE = blkdiagExtend(W_RF_PE);
        W_Deq_PE = blkdiagExtend(W_D_PE);
        desymbolPE = beta_PE^(-1)*W_Deq_PE'*W_RFeq_PE'*symbol_rPE;
        desymbolPE = desymbolPE(:);
        symbol_rPE = step(hDemod,desymbolPE);
        biterrPE(n) = length(find(data~=symbol_rPE));
%         if biterr_FD(n)>biterrPE(n)
%             save A.mat
%             pause;
%         end
        
        if mod(n,100)==0
            fprintf('at %d dB, iteration : %d times  ',SNR_dB(snr_idx), n ) ;
            disp(datestr(now)); 
        end

    end
    BitErr(snr_idx) = sum(biterr)./(N_sym*N_s*N_u*2*N_loop)
    BitErr_FD(snr_idx) = sum(biterr_FD)./(N_sym*N_s*N_u*2*N_loop)
%     BitErrY(snr_idx) = sum(biterrY)./(N_sym*N_s*N_u*2*N_loop)
    BitErrPE(snr_idx) = sum(biterrPE)./(N_sym*N_s*N_u*2*N_loop)
    
    sumMSE(snr_idx) = sum(mse)./(N_loop)
    sumMSE_FD(snr_idx) = sum(mse_FD)./(N_loop)
%     sumMSEY(snr_idx) = sum(mseY)./(N_loop)
    sumMSEPE(snr_idx) = sum(msePE)./(N_loop)
end
figure
semilogy(SNR_dB,BitErr_FD,'--');
hold on;
semilogy(SNR_dB,BitErr,'-*b');
% semilogy(SNR_dB,BitErrY,'-+')
semilogy(SNR_dB,BitErrPE,'-<');
legend('Fully Digital','Proposed Scheme','PE-AltMin');

% title(['Downlink'] );
xlabel('SNR(dB)'); ylabel('BER');
grid on;

figure
semilogy(SNR_dB,sumMSE,'-*b');
hold on;
semilogy(SNR_dB,sumMSE_FD,'--');
% semilogy(SNR_dB,sumMSEY,'-+')
semilogy(SNR_dB,sumMSEPE,'-<');
legend('Uni','Digital','PE');
title(['Downlink'] );
xlabel('SNR(dB)'); ylabel('MSE');
grid on;
disp(datestr(now));                % show time