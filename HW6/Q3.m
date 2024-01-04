clear all, close all, clc

u = @(x,y) exp(x*y) + x^2; 

x = 1;          
y = 2;
hvals = 2.^-(2:12);

MY_deltaU = FivePt_TwoD_Laplacian(u, x, y, hvals); 

Df_X =  (y^2 * exp(x*y)) + 2;
Df_Y =  x^2 * exp(x*y);

EXACT_deltaU =  Df_X + Df_Y;

for i = 1:length(hvals)                                              
    Rel_Enorm(i) = norm(MY_deltaU(i) - EXACT_deltaU) / norm(EXACT_deltaU);
end

% LOG-LOG plot of ||e|| VS. h 
loglog(hvals, Rel_Enorm, 'b-o')            
hold on
C = 1/24;                           
y = C * hvals(3:5) .^ 2;                 
loglog(hvals(3:5), y, 'r-o')
xlabel('h')
ylabel('Relative Error Norm')
title('Log-Log plot of ||e|| vs. h ')
legend('Relative Error Norm', '\Delta h^2')
set(legend, 'Location', 'NorthWest')



% Determine the slope of the TINY SEGMENT for hvals(3:5) 
y = @(x) C * x .^ 2;

y2 = log(y(0.0625));
y1 = log(y(0.0156));
  
x2 = log(hvals(3));
x1 = log(hvals(5));


slope = (y2 - y1) / (x2 - x1)

function [deltaU] = FivePt_TwoD_Laplacian(u, x, y, hvals) 
    for i = 1:length(hvals)  
        deltaU(i) = (u(x-hvals(i), y) + u(x+hvals(i), y) + u(x, y-hvals(i)) + u(x, y+hvals(i)) - 4*u(x, y)) / hvals(i)^2;
    end
end