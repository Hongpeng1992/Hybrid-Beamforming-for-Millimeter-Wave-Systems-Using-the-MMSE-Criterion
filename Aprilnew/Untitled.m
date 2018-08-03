flag = 0
n = 10;
for x = 1:10000
A =rand(n,n);
A = A*A';
B = rand(n,n);
B = B*B';
a=eig(A);
c = eig(A+B);
for i = 1:n
    if c(i)<0
        flag = 1
    end
end
end