function [cost, storedb] = getCost1(problem, x, storedb,H,Vn,N_RF)
                [~,N_s,N_k]=size(H); 
                N_r = size(x,1)/N_RF;
                x = reshape(x,N_r,N_RF);
                cost = zeros(1,N_k);
                for i = 1:N_k
                cost(i) = trace((H(:,:,i)'*x*(x'*x)^(-1)*x'*H(:,:,i)/Vn(i)+eye(N_s))^(-1));
                end
                cost = sum(cost);
            
end