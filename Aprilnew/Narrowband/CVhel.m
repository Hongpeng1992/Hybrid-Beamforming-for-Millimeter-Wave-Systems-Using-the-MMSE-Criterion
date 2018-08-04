function [X] = CVhel(X)
[a,b]=size(X);
for i = 1:a
    for m = 1:b
        if X(i,m) == 0 
            X(i,m) = X(i,m-1);
        end
    end
end