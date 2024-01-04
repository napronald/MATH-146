clear all; clc; close all

A = [4 0 0; -2 6 1; -1 1 7];
b = [3;9;-6];

tol= 1e-12;
max_iter = 10;

% A = [4 0 0; -2 6 1; -1 1 7];
A = [3 -2 1; 1 3 2; -1 2 4]
b = [3;9;-6];
A\b
% tolerence = 1e-12;
% MaxNumOfIter = 25;

% x=GaussSeidel(A,B,tolerence,MaxNumOfIter)

[x,iter] = GaussSeidel(A,b,tol,max_iter)

function [x,iter] = GaussSeidel(A,b,tol,max_iter)
    [n,m]=size(A);
    
    iter=1; % actual number of iterations used
    x0=zeros(n,1); % we need a zero matrix to build in
    x=zeros(n,1); % initialise output matrix x
    
    while (iter<max_iter) % loop ends when exceed max no. of iterations 
        for i=1:n
            I=[1:i-1 i+1:n];
            x(i)=(b(i)-A(i,I)*x(I))/A(i,i);
        end
        esp=(abs(x(i)-x0)/abs(x(i)));
        if max(esp)<tol
           break;
        end
        x0=x; % assign the new x to x0
        iter=iter+1;% no. of iterations
    end
end

% function x=GaussSeidel(A,B,tolerence,MaxNumOfIter)
%     % A and B are input matrix, err is the tolerence, NumOfIter is the
%     % number of iterations
%     % x is the output solution
%     % take the size of two matrix
%     [n,m]=size(A);
%     [u,~]=size(B);
%     
%    % check if A is a square matrix and weather dimensions of A and B match or not
%     if n~=m
%         error('Matrix A must be a square matrix');
%     elseif  n~=u
%         error('The number of rows of A must be the same as that of B');
%     end
%     
%     counter=1; % actual number of iterations used
%     x0=zeros(n,1); % we need a zero matrix to build in
%     x=zeros(n,1); % initialise output matrix x
%     
%    %check for conditions of GaussSeidel, see if the matrix is strictly
%    %diagonally dominant or not
%     for i=1:n
%     s=0;
%     %sum of the row
%     for j=1:n
%         if i~=j 
%             s=s+abs(A(i,j));
%         end
%     end
%     
%     % invalid when A(1,1)(2,2)...(n,n)>sum of the related row
%     % OR any diagnal entry is zero which means A is not a strictly
%     %diagonally dominant matrix or a positive definite matrix
%     if s<abs(A(i,j))
%         fprintf('The conditions of Guass-Seidel have not met\nPlease use other methods to find solution matrix x\n');
%         return;
%     end
%     end
%        
%     
%     while (counter<MaxNumOfIter) % loop ends when exceed max no. of iterations 
%         %Gauss-Seidel Iteration
%         for i=1:n
%             I=[1:i-1 i+1:n];
%             x(i)=(B(i)-A(i,I)*x(I))/A(i,i);
%         end
%         %calculate error and compare with tolerence entered
%         esp=(abs(x(i)-x0)/abs(x(i)));
%         %break if good enough
%         if max(esp)<tolerence
%            break;
%         end
%         x0=x; % assign the new x to x0
%         counter=counter+1;% no. of iterations
%     end
%     disp(['The number of iterations is ', num2str(counter)]);
%     if counter>=MaxNumOfIter
%     fprintf('Maximum number of iterations has exceeded\n');
%     end
% end
