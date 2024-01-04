clear all
clc

n = 3;
A = 20*rand(n,n);

A = [1,2,3;4,5,6;7,8,0]

AUG = [A,eye(size(A))]
for j=1:length(A)
    AUG(j,:) = AUG(j,:) / AUG(j,j)
    for i = 1:length(A)
        if i ~=j
            disp('runs')
            AUG(i,:) = AUG(i,:) - AUG(i,j) * AUG(j,:)
        end
    end
%     inv = AUG(:,size(A)+1:size(A)*2);
end
