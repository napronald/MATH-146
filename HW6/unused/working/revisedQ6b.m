clear all; clc;

n=1000;
[A,b] = nap.rand_sdd_matrix(n);

tolerence = 1e-12;

MaxNumOfIter = 25;

[x,iter] = GaussSeidel(A,b,tolerence,MaxNumOfIter,n);
x_matlab = A\b;

iter
relative_difference = norm(x-x_matlab) / norm(x)

function [x,iter] = GaussSeidel(A,b,tol,MaxNumOfIter,n)
    iter=1; 
    x0=zeros(n,1); 
    x=zeros(n,1); 
    
    while (iter<MaxNumOfIter) 
        for i=1:n
            I=[1:i-1 i+1:n];
            x(i)=(b(i)-A(i,I)*x(I))/A(i,i); 
        end
        
        error = abs(x - x0);
        if norm(error)<tol
           break;
        end
        x0=x; 
        iter=iter+1;
    end
end