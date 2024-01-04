clear all
clc
close all

% Plot of step size=1 provided in files

step = 1; % Adjust step size
% Note:step=1 will take around 1 minute to run since our algo is O(n^2) 
% with increasingly larger matrices

for n = 5:step:400
    A = 20*rand(n,n); % Formulating random A,z
    z = 20*rand(n,1);
    
    B = A*z;
    [x] = nap.Ax_B(A,B); % using function from nap.m file
    
    rel_1(n) = norm(x-z,2) / norm(z,2); % Computing rel_err + cond_num
    cond_num_1(n) = cond(A);
end

% plotting first batch of points
figure(1)
loglog(rel_1,cond_num_1,'.')
grid on 
hold on

for n = 401:step:500
    A = 20*rand(n,n); % Formulating random A,z
    z = 20*rand(n,1);
    
    B = A*z; 
    [x] = nap.Ax_B(A,B); % using function from nap.m file
    
    rel_2(n) = norm(x-z,2) / norm(z,2); % Computing rel_err + cond_num
    cond_num_2(n) = cond(A);
end

% plotting second batch of points
loglog(rel_2,cond_num_2,'r.')
xlabel('Relative Error')
ylabel('Condition Number')
title('Relative Error vs Condition Number')
legend('n = 5 to 400','n = 401 to 500')