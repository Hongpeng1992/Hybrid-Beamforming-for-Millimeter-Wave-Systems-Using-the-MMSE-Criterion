function [ W_D,W_RF ] = WEVD( H, N_RF ,Vn )

[Nt, Ns,N_k] = size(H);
X = zeros(Nt,Nt,N_k);
W_D = zeros(N_RF,Ns,N_k);
t1 = clock;
for k = 1:N_k
    X(:,:,k) = ( 1/Vn(k)/Nt*H(:,:,k)*H(:,:,k)');
end
M =sum(X,3);
t2 = clock;
a = etime(t2,t1);
t1 = clock;
for i = 1:N_RF/2
[V,D] = eig(M);
[~,IX] = sort(diag(D),'descend');

    eig_Vector = V(:,IX);
    W_RF = eig_Vector(:,1:N_RF);
    W_RF = exp(1i*angle(W_RF));
    t2 = clock;
    b = etime(t2,t1);
end
t2 = clock;
t = etime(t2,t1);
t1 = clock;
for i = 1:N_k
    M = W_RF'*(H(:,:,i))*H(:,:,i)'*W_RF+ Vn(i)*(W_RF)'*W_RF;
    W_D(:,:,i) = inv(M)*W_RF'*H(:,:,i);
end
t2 = clock;
c = etime(t2,t1);