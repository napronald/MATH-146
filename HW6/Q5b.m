clear all; clc; close all

f = @(x,y) -exp((x-0.25).^2 + (y-0.25).^2); % Given f

uold = 0; % Initialize uold
for i=2:10 
N = (2^i)-1; % for N = 2^3 to 2^10
h = 1 / (N + 1); 

% Interior Gridpoints
xg=h*(1:N);
yg=h*(1:N); 
[xg,yg]=ndgrid(xg,yg);

%Compute right hand side
bmat = f(xg,yg); 
b = reshape(bmat,N*N,1);

% Compute A matrix
[L2] = nap.lap2d(N,N);
A = (1/h^2)*L2;

% Solve for u and format into vector
unew = A\b;
unew = reshape(unew,N,N);

% skip first iteration, no other point to compare with
if i > 1
    unew_skip = unew(2:2:end,2:2:end); % skip grid points

    diff = unew_skip - uold; % Compute difference
    diffmax(i-1) = max(abs(diff(:))); % Compute Max Norm
end
uold = unew; % save u for next iteration

end

figure(1)
mesh(xg,yg,unew) % 3D plot

figure(2)
N = (2.^(2:10))-1; % Given N's

loglog(N(2:end),diffmax(2:9),'-o') % plot excluding index 1
hold on

% plot the slope of the TINY SEGMENT for N(4:6) 
C = 1/24;                                 
y = C ./ (N(4:6) .^ 2);             
loglog(N(4:6), y, 'r-o') 

xlabel('N')
ylabel('Max Norn')
title('Log-Log plot of ||e|| vs. N')
legend('Max Norm', '\Delta h^2')
set(legend, 'Location', 'NorthEast')

% Determine the slope of the TINY SEGMENT for hvals(3:5) 
y = @(x) C / (x .^ 2)  ;

y2 = log(y(64))        ; 
y1 = log(y(256))       ; 

x2 = log(N(4))     ; 
x1 = log(N(6))     ; 

slope = (y2 - y1) / (x2 - x1)