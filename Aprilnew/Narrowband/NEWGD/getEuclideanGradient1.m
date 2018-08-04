function [egrad, storedb] = getEuclideanGradient1(problem, x, storedb,H,Vn,N_RF)

                store = getStore(problem, x, storedb);
                if ~isfield(store, 'egrad__')
                    store.egrad__ = newGD(H,Vn,N_RF,x);
                    storedb = setStore(problem, x, storedb, store);
                end
                egrad = store.egrad__;
         
        end

   
