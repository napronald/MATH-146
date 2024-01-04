function [x] = back_sub(U,B)
n = length(B);
x = zeros(length(B),1);

x(n) = B(n)/U(n,n); 
for i = n-1:-1:1 
    sum = 0;
    for j = n:-1:1
        sum = sum + U(i,j) * x(j);
    end
    x(i) = (B(i) - sum) / U(i,i);
end
end
