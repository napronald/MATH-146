clear all; clc; close all


u_actual = @(x) x - sin(pi*x) / pi^2;               % ACTUAL solution


% Given boundary conditions (B.C.s: u(0) = 0 ----> x_0 = a  &  u(1) = 1 ----> x_N = b) 
a = 0;
b = 1;


% For values of N = 2^3 TO 2^12 
Nvals = 2.^(3:12);

for j=1:10
h = (b - a) / (Nvals(j) + 1);           % Nvals = 2^3 TO 2^12 =====> h = 10 terms TOTAL 


% For [a,b] ----> from x_0 = a TO x_N = b ====> x = a TO b w/ increments of h 
x = a:h:b;


% Define matrix A 
A = (1/h^2) * (-2*eye(Nvals(j)) + diag(ones(1, Nvals(j)-1),  1) + diag(ones(1, Nvals(j)-1),  -1));                                                                                               
    % Components of Vector a ---> from 2 TO n  OR  from 1 TO n-1   =====> Vector a: 1's on SUB-Diagonal   ---> desginated @ -1  (shifts DOWN from the diagonal)
    % Components of Vector b ---> from 1 TO n                      =====> Vector b: -2's on DIAGONAL      ---> desginated @ 0                                                                                     
    % Components of Vector c ---> from 1 TO n-1                    =====> Vector c: 1's on SUPER-Diagonal ---> desginated @ 1   (shifts UP from the diagonal) 


% Solve for B
f = sin(pi*x(2:end-1))';

B(1) = f(1)- (a/h^2);                       % 1st component of Vector B =====> f_1 - a/h^2 
for i = 2:length(f)-1
    B(i) = f(i);
end
B(length(f)) = f(length(f)) - (b/h^2);      % LAST component of Vector B ====> f_N - b/h^2 


% Solve Au = B ----> A : Matrix  &  B : Vector 
u = A \ B'  ;                                         


% Check this part 
for k = 1:length(u)
    approx(k+1) = u(k);
end
approx(1) = a;
approx(end+1) = b;
approx;

Rel_Enorm(j) = norm(approx-u_actual(x), inf) / norm(u_actual(x), inf);            % Relative Error Norm w.r.t MAX-NORM 
end

% Plot the APPROXIMATED solutions & the TRUE solution 
subplot(2,1,1)
plot(x, u_actual(x), 'bo-', x, approx, 'rx-', 'MarkerSize', 5);
xlabel('x');
ylabel('f(x)')
title('Comparison between the Approximated solutions and the True solution')
legend('True solution', 'Approximated solutions')
set(legend, 'Location', 'NorthWest')




% LOG-LOG plot of ||e|| VS. N 
subplot(2,1,2)
loglog(Nvals, Rel_Enorm, 'b-o')            
hold on
C = 1/24;                                       % Nvals = C / N^2 
y = C ./ (Nvals(4:6) .^ 2);                     % P = 2 ---> Expect 2nd-order accuracy 
loglog(Nvals(4:6), y, 'r-o')                         % Plots the TINY SEGMENT 
xlabel('N')
ylabel('Relative Error Norm')
title('Log-Log plot of ||e|| vs. N')
legend('Relative Error Norm', '\Delta h^2')
set(legend, 'Location', 'NorthEast')


% Determine the slope of the TINY SEGMENT 
y = @(x) C / (x .^ 2)  ;

y2 = log(y(64))        ; 
y1 = log(y(256))       ; 

x2 = log(Nvals(4))     ; 
x1 = log(Nvals(6))     ; 


slope = (y2 - y1) / (x2 - x1)