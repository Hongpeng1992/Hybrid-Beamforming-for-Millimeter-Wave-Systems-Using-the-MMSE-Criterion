function [ H,W_codebook,F_codebook ] = ChannelULA( N_r,N_t,N_u,L )
%  ULA信道产生
%  对比曲线所需码本
H = zeros(N_r, N_t,N_u);                 %信道
H_path = zeros(N_r, N_t,L);
W_codebook = zeros(N_r,L,N_u);         %用户端RF码本
F_codebook = zeros(N_t,L,N_u);         %基站端RF码本
for r = 1:N_r
    %-----------AOA MS
    aoa = 2*pi* rand(N_u,L);          %服从（0,2*pi）的均匀分布
    aoa = sin(aoa);
    %-----------AOD BS
    aod = 2*pi* rand(N_u, L);         %服从（0,2*pi）的均匀分布
    aod = sin(aod);
    %-----------Complex path gain
    alpha = complex(randn(N_u,L),randn(N_u,L))/sqrt(2);
    
    signature_t = [0:(N_t-1)]';
    signature_t = 1i*pi* signature_t;               %为接下来的signature构造做准备
    signature_r = [0:(N_r-1)]';
    signature_r = 1i*pi* signature_r;               %为接下来的signature构造做准备

    for K_i = 1:N_u
        for L_i= 1: L
            H_path(:,:,L_i)=alpha(K_i,L_i)*exp(aoa(K_i,L_i)*signature_r)*exp(aod(K_i,L_i)*signature_t)';
            W_codebook(:,L_i,K_i) = exp(aoa(K_i,L_i)*signature_r);
            F_codebook(:,L_i,K_i) = exp(aod(K_i,L_i)*signature_t);
        end
        H(:,:,K_i) = sqrt(N_t*N_r/L)*sum(H_path,3);
    end
end

end

