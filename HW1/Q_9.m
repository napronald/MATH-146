clear all
clc
close all

A = 30.*rand(1000,40);
B = 30.*rand(40,4000);

disp('self built function')
[M] = mat_mult(A,B);

disp('Built in function')
tic
V = A*B;
toc

% validation = isequal(fix(V),fix(M))

function [M] = mat_mult(A,B)
M=zeros(size(A,1),size(B,2));
tic
if size(A,2) == size(B,1) 
    for i=1:size(A,1)
        for j=1:size(B,2)
            for k=1:size(A,2)
                M(i,j) = M(i,j) + A(i,k) * B(k,j);
            end
        end
    end
    M;
else
    disp('Error, matrices cannot be multiplied.')
end
toc
end