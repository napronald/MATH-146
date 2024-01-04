clear all, close all, clc

% INPUTS: 
u = @(x,y) exp(x*y) + x^2  ;            % Given function 

x = 1  ;          
y = 2  ;
hvals = 2.^-(2:12);                     % Values of h: 2^-2 TO 2^12 

MY_deltaU = FivePt_TwoD_Laplacian(u, x, y, hvals)  ;        % APPROXIMATE solution, using my function 

Df_X =  (y^2 * exp(x*y)) + 2   ;              % 2nd Partial Derivative of the function (u) w.r.t X
Df_Y =  x^2 * exp(x*y)         ;              % 2nd Partial Derivative of the function (u) w.r.t Y

EXACT_deltaU =  Df_X + Df_Y    ;              % ACTUAL solution = y^2exp(xy) + x^2exp(xy) + 2 

% VALIDATE CODE
for i = 1:length(hvals)                                              % Relative Errors for values of h (2^-2 TO 2^-12) 
    Rel_Enorm(i) = norm(MY_deltaU(i) - EXACT_deltaU) / norm(EXACT_deltaU)  ;
end

x = hvals(1:end-1)   ;
y = Rel_Enorm(2:end) ;
logx = log(x)  ;            % Take the log of x (or h) 
logy = log(y)  ;            % Take the log of y (or error) 

% Plot a small segment of the solution ----> Determine the
x1 = logx(3)  ;             % end-1: 2nd to last point 
y1 = logy(3)  ;
x2 = logx(4)    ;             % end: Last point 
y2 = logy(4)    ;

xpts = [x1-1 x2-1]  ;
ypts = [y1-1 y2-1]  ;

% Determine the SLOPE of the SMALL SEGMENT in the LOG-LOG plot 
slope =  (y2 - y1) / (x2 - x1)

% LOG-LOG plot of ||e|| VS. h 
% loglog(hvals, Rel_Enorm, 'bo')

% Plot the LOG-LOG plot of the RELATIVE ERRORS 

hold on

loglog(hvals, Rel_Enorm, 'r-o')
% hold on
% plot(xpts,ypts,'b-o')
% xlabel('h')
% ylabel('Relative Error Norm')
% title('Log-Log plot of ||e|| vs. h ')
% legend('Relative Error Norm', '\Delta h^2')



%% Function to use 5-POINT 2-D LAPLACIAN  
function [deltaU] = FivePt_TwoD_Laplacian(u, x, y, hvals) 
    % INPUTS: % u ---> FUNCTION
              % x
              % y
              % hvals 

    for i = 1:length(hvals)  
        deltaU(i) = (u(x-hvals(i), y) + u(x+hvals(i), y) + u(x, y-hvals(i)) + u(x, y+hvals(i)) - 4*u(x, y)) / hvals(i)^2;
    end
end
