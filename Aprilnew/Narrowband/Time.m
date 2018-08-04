clear all; close all; clc
disp(datestr(now));                % show time
%---------structure parameter---------%
SNR_dB = [-16];
SNR = 10 .^ ( SNR_dB/10 );
P = 1;
N_t = 64;                  % Number of transmit antennas
N_r = 64;                  % Number of receive antennas
N_RF = 4;                  % Number of RF chain
Ns = 2;                    % Number of data streams
N_s = 64;                  % Number of symbols for one channel realization
fprintf('simulate the %d * %d * %d performance of BER at SNR = %d : %d \n ', N_t,N_r,N_RF, SNR_dB(1) , SNR_dB(end) );
%%---channel parameter----------%
N_loop =500;          % Number of loop per SNR
cn = 3;                    %comparing varieties
%%---data matrices----------%
BitErr = zeros(cn+1,length(SNR_dB));     %BER performance
Jiankong = zeros(1,length(SNR_dB));   % optimal works in partial channel
Diedaicishu = zeros(2,length(SNR_dB));
% load('H.mat')
% load('AT.mat')
% load('AR.mat')
for snr_idx = 1: length(SNR_dB)
    %%---data matrices----------%
    stanMSE = zeros(cn+1 ,N_loop);            %record the MSE
    biterr = zeros(cn+1 ,N_loop);
    dd1 = 0;
    dd2 = 0;
    Vn = 1 / 10^(SNR_dB(snr_idx)/10);      % Noise Power
    t1 = clock;
    n = 1;
    SMSE1 = zeros(N_loop,11);
    SMSE3 = zeros(N_loop,11);
    SMSE2 = zeros(N_loop,11);
    Optime = 0;
    OMPtime = 0;
    GEtime = 0;
    MOtime = 0;
    yuweitime = 0;
    %MSE1 = zeros(N_loop,20);
    %MSE = zeros(N_loop,20);
    for n = 1 : N_loop
        %% 先产生信
        [H,AT,AR]  = OMPH(N_t,N_r);
        %          H = h(:,:,n);
        %          AT = at(:,:,n);
        %          AR = ar(:,:,n);
        %% obtain beamforming matrix
        
        [V_opt, W_opt, MSE_opt] = MSEopt( H, Vn ,Ns);     % the optimal
        if ( norm(V_opt(:,end))==0 )
            Jiankong(snr_idx) =    Jiankong(snr_idx) + 1;              %avoid the bad case
            continue;
        end
        t1 = clock;
        [~] =  MSEopt( H, Vn ,Ns);
        t2 = clock;
        Optime = Optime +etime(t2,t1);
        t1 = clock;
        [T(1,n)] =TGEini(Ns,N_RF,H,Vn,W_opt);
        t2 = clock;
        GEtime = GEtime + etime(t2,t1);
        t1 = clock;
        %[T(2,n)] = TMOMSEini(Ns,N_RF,H,Vn,W_opt);
        t2 = clock;
        MOtime = MOtime +etime(t2,t1);
        t1 = clock;
        [T(3,n)] = TOMPini(Ns,N_RF,H,Vn,W_opt,AT,AR);
        t2 = clock;
        OMPtime = OMPtime +etime(t2,t1);
        t1 = clock;
        [T(4,n)]=TYUWEI(Ns,N_RF,H,Vn);
        t2 = clock;
        yuweitime = yuweitime +etime(t2,t1);
        %% simulate the transmission
        
        %% display
        if (n==10)
            mytoc(N_loop,t1);
        end
    end
    BitErr(:,snr_idx) = sum(biterr,2)./(N_s*Ns*2*(N_loop-Jiankong(snr_idx)))
    getime = GEtime/N_loop
    motime = MOtime/N_loop
    optime = Optime/N_loop
    omptime = OMPtime/N_loop
    X = sum(T,2)/N_loop
end
disp(datestr(now));
