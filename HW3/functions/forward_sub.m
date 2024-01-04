function [x] = forward_sub(L,B)
n = length(B);
x = zeros(length(B),1);

x(1) = B(1)/L(1,1); 
for i = 2:n 
    sum = 0;
    for j = 1:i-1 
        sum = sum + L(i,j) * x(j); 
    end
    x(i) = (B(i) - sum) / L(i,i);
end
end

