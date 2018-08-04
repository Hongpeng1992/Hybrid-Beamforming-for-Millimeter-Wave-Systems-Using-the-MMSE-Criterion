function [t] = TYUWEI(Ns,N_RF,H,Vn)
[Nr,Nt,Nk]  = size(H);
V_RF = yuweiA1(Vn,H,N_RF,Nt);
Q = (V_RF'*V_RF);
T = Q^(-0.5);
for k = 1:Nk
    L = H(:,:,k)*V_RF*T;
    [~,D,V] = svd(L);
    [~,IX] = sort(diag(D),'descend');
    M = V(:,IX);
    U = M(:,1:Ns);
    V_D(:,:,k) = T*U;
    V_D(:,:,k) = V_D(:,:,k)/norm(V_RF*V_D(:,:,k),'fro');
end
t1 = clock;
W_RF  = yuweiA2(Vn,H,N_RF,Nr,V_D,V_RF);
t2 = clock;
t = etime(t2,t1);
for k = 1:Nk
    J = W_RF'*H(:,:,k)*V_RF*V_D(:,:,k)*V_D(:,:,k)'*V_RF'*H(:,:,k)'*W_RF+Vn*W_RF'*W_RF;
    W_D(:,:,k) = J^(-1)*W_RF'*H(:,:,k)*V_RF*V_D(:,:,k);
end
    