clear all
clc

A = [1,2,3,1;4,5,6,2;7,8,0,4;0,1,3,1];

[L, U, P] = LU_pivot(A)

function [L, U, P]=LU_pivot(A)

n=length(A);
L=eye(n); 
P= eye(n);
U=A;

for j=1:n-1
    [max_num, index] = max(abs(U(j:n,j)));
    index=index+j-1
    if index~=j
        temp=U(j,:);
        U(j,:)=U(index,:);
        U(index,:)=temp;
        temp=P(j,:);
        P(j,:)=P(index,:);
        P(index,:)=temp;
        if j >= 2
            temp=L(j,1:j-1);
            L(j,1:j-1)=L(index,1:j-1);
            L(index,1:j-1)=temp;
        end
    end

    for i=j+1:n
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;
        U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
    end
end
end

%%
clear all
clc

A = [1,2,3,1;4,5,6,2;7,8,0,4;0,1,3,1];

[L, U, P] = LU_pivot(A)

function [L, U, P]=LU_pivot(A)

n=length(A);
L=eye(n); 
P= eye(n);
U=A;

for j=1:n-1
    [max_num, index] = max(abs(U(j:n,j)));
    index = index+j-1;
    if index ~= j
        U([index, j], :) = U([j, index], :);
        P([index, j], :) = P([j, index], :);

        if j >= 2
%             temp=L(j,1:j-1)
%             L(j,1:j-1)=L(index,1:j-1)
%             L(index,1:j-1)=temp;

%             L([1:j-1, index], :) = L([index, 1:j-1], :);
%             L([j, 1:j-1], :) = L([index, 1:j-1], :)
            L([index, j], :) = L([j, index], :)

        end
    end

    for i=j+1:n
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;;
        U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
    end
end
end


%%
clear all
clc
close all

rel_1 = zeros(1,485);
cond_norm_1 = zeros(1,485);

for n = 5:490
    A = 20*rand(n,n);
    z = 20*rand(n,1);
    
    B = A*z;

    [x] = nap.Ax_B(A,B);
    
    rel_1(n-4) = norm(x-z,2) / norm(z,2);
    cond_num_1(n-4) = cond(A);

end


figure(1)
loglog(rel_1,cond_num_1,'.')
grid on 
hold on

rel_2 = zeros(1,10);
cond_norm_2 = zeros(1,10);

for n = 491:500
    A = 20*rand(n,n);
    z = 20*rand(n,1);
    
    B = A*z; 

    [x] = nap.Ax_B(A,B);
    
    rel_2(n-490) = norm(x-z,2) / norm(z,2)
    cond_num_2(n-490) = cond(A)

end

loglog(rel_2,cond_num_2,'r.')
xlabel('Norm')
ylabel('condition number')
title('Norm vs Condition Number')
legend('5 to 490','490 to 500')
















