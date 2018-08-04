function [ W_RF,V_RF,W_D,V_D,beta,mse2 ] = Uni_GED_Opt( N_r,N_t,N_RFt,N_RFr,N_u,N_s,Vn,P,H ,W_opt)
%  预编码矩阵优化
%  基于广义特征值分解
%  发送端和接收端迭代优化
%  接收端数字阵为酉阵

V_RF = randn(N_t,N_RFt)+1i*randn(N_t,N_RFt);
phase = angle(V_RF);
V_RF = exp(1i*phase);
W_RF = W_opt;
% W_RF = randn(N_r,N_RFr,N_u)+1i*randn(N_r,N_RFr,N_u);  %随机初始点
% phase = angle(W_RF);
% W_RF = exp(1i*phase);
W_RFeq = blkdiagExtend(W_RF);
W_D = zeros(N_RFr,N_s,N_u);
for k_i = 1:N_u
    W_D(:,:,k_i) = eye(N_s);
end
W_Deq = blkdiagExtend(W_D);
Weq = W_RFeq*W_Deq;
Heq = ChannelEqu( H );
V_D = inv((V_RF'*Heq'*W_RFeq*W_Deq*W_Deq'*W_RFeq'*Heq*V_RF+Vn*trace(W_Deq'*W_RFeq'*W_RFeq*W_Deq)/P*V_RF'*V_RF))*V_RF'*Heq'*W_RFeq*W_Deq;
beta = (trace(V_RF*V_D*V_D'*V_RF')/P)^(-1/2);
mse2 = norm(eye(N_u*N_s)-W_Deq'*W_RFeq'*Heq*V_RF*V_D,'fro')+...
    beta^(-2)*Vn*trace(W_Deq'*W_RFeq'*W_RFeq*W_Deq);
mse = 0;
% for N_loop = 1:10
while abs(mse-mse2)>10^(-2)
    mse = mse2;

    % Optimazing V_RF based on MMSE
    for m = 1:5            % N_loop=5 convergence
        for j = 1:N_RFt
            V_RF_tilde = V_RF;
            V_RF_tilde(:,j) =[];               %去掉第j列的V_RF
            Aj = (1/(Vn*N_t*N_r*N_u*N_s/P))*W_Deq'*W_RFeq'*Heq*V_RF_tilde*V_RF_tilde'*Heq'*W_RFeq*W_Deq+eye(N_u*N_RFr);
            A = (1/(Vn*N_t*N_r*N_u*N_s/P))*Heq'*W_RFeq*W_Deq*(Aj^(-2))*W_Deq'*W_RFeq'*Heq;
            B = 1/N_t*eye(N_t)+(1/(Vn*N_t*N_r*N_u*N_s/P))*Heq'*W_RFeq*W_Deq*(Aj^(-1))*W_Deq'*W_RFeq'*Heq;
            [V,D] = eig(A,B);
            [newD,IX] = sort(diag(D),'descend'); %特征值降序排列
            eig_Vector = V(:,IX);                  %取出对应特征值的特征向量
%             V_RF(:,j) = eig_Vector(:,1);
            phase = angle( eig_Vector(:,1) );
            V_RFjopt = exp(1i*phase);            %取相位，保证模为1
            V_RF(:,j) = V_RFjopt;
        end
    end
    % Optimize V_D
    V_D = inv((V_RF'*Heq'*W_RFeq*W_Deq*W_Deq'*W_RFeq'*Heq*V_RF+Vn*trace(W_Deq'*W_RFeq'*W_RFeq*W_Deq)/P*V_RF'*V_RF))*V_RF'*Heq'*W_RFeq*W_Deq;
    beta = (trace(V_RF*V_D*V_D'*V_RF')/P)^(-1/2);
    
    % Optimazing W_RF based on MMSE
    for m = 1:5            % N_loop=5 convergence
        for Ki = 1:N_u
            for NRFi = 1:N_RFr
                j = (Ki-1)*N_RFr+NRFi;
            W_RF_tilde = W_RFeq;
            W_RF_tilde(:,j) =[];               %去掉第j列的W_RF
            Aj = V_RF'*Heq'*W_RF_tilde*W_RF_tilde'*Heq*V_RF+(Vn*N_t*N_r*N_u*N_s/P)*eye(N_RFt);
            A = H(:,:,Ki)*V_RF*(Aj^(-2))*V_RF'*H(:,:,Ki)';
            B = 1/N_r*eye(N_r)+H(:,:,Ki)*V_RF*(Aj^(-1))*V_RF'*H(:,:,Ki)';
            [V,D] = eig(A,B);
            [newD,IX] = sort(diag(D),'descend'); %特征值降序排列
            eig_Vector = V(:,IX);                  %取出对应特征值的特征向量
%             W_RF(:,:,j) = eig_Vector(:,1);
            phase = angle( eig_Vector(:,1) );
            W_RFjopt = exp(1i*phase);            %取相位，保证模为1
            W_RF(:,NRFi,Ki) = W_RFjopt;
            W_RFeq = blkdiagExtend(W_RF);
            end
        end
    end
    % Optimize W_D
    W_D = WDopt( N_RFr,N_s,N_u,W_RF,H,V_RF,V_D );
    W_Deq = blkdiagExtend(W_D);
    mse2 = norm(eye(N_u*N_s)-W_Deq'*W_RFeq'*Heq*V_RF*V_D,'fro')+...
            beta^(-2)*Vn*trace(W_Deq'*W_RFeq'*W_RFeq*W_Deq);
end
V_D = beta*V_D;
end

