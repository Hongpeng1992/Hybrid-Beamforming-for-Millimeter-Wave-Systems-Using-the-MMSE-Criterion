function [H] = Channel(N_t,N_r)
L=3;
E_aoa = 2*pi* rand(L,1);                               %服从（0,2*pi）的均匀分布
aoa = cos(E_aoa);
%-----------AOD
E_aod = 2*pi* rand(L, 1);                               %服从（0,2*pi）的均匀分布
aod = cos(E_aod);
signature_t = [0:(N_t-1)]';
signature_t = 1i*pi* signature_t;                           %为接下来的signature构造做准备
signature_r = [0:(N_r-1)]';
signature_r = 1i*pi* signature_r;                           %为接下来的signature构造做准备

H_ray = zeros(N_r, N_t,L);

for m = 1: L
    H_ray(:,:,m)=complex(randn(1),randn(1))/sqrt(2)*exp((aoa(m)*signature_r))*exp((aod(m)*signature_t))'/sqrt(N_t*N_r);
end

H = sqrt(N_t*N_r/L)*sum(H_ray,3);
