function [cost, storedb] = getCost1(problem, x, storedb,H,Vn,N_RF)
                [~,N_s]=size(H); 
                N_r = size(x,1)/N_RF;
                x = reshape(x,N_r,N_RF);
                cost = trace((H'*x*(x'*x)^(-1)*x'*H/Vn+eye(N_s))^(-1));
              
            
end
