function [diedai] = march4(X)
N = size(X,1);
X = X(:);
count = 0;
for i = 1:length(X)
    if X(i)>0
        count = count+1;
    end
end
diedai = count/N;