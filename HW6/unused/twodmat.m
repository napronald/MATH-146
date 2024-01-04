clear all; clc; close all

% N = 3
% M = 4
% 
% dx = 1/N
% dy = 1/M
% 
% A = zeros((M-1)*(N-1));
% A = A + eye((M-1)*(M-1)) * (-2/dx^2-2/dy^2);
% 
% for i=1:N-1
%     A((i-1)*(M-1)+1:i*(M-1),(i-1)*(M-1)+1:i*(M-1)) = A((i-1)*(M-1)+1:i*(M-1),(i-1)*(M-1)+1:i*(M-1)) + diag(ones(M-2,1) / dy^2,1)
% end


u = @(x,y) -2*(x-y)*exp((y-0.25)^2 - (x-0.25)^2);

N = 3;
a=0;
b=1;

h = 1 / (N + 1);
x = a:h:b;
y = a:h:b;

xg=h*(1:N); 
yg=h*(1:N); 

[xg,yg]=ndgrid(xg,yg);

x_int = xg + h;
y_int = yg + h;

A = zeros(N*N);
f = zeros(N*N,1);


% Loop over interior grid points
for i=1:N
    for j=1:N
        % Index of current grid point
        ind = (j-1)*N + i;
        
        % Finite difference approximation of Laplacian
        A(ind,ind) = -4/h^2;
        if i > 1
            A(ind,ind-1) = 1/h^2;
        else
            % Boundary condition
            f(ind) = f(ind) - 0/h^2;
        end
        if i < N
            A(ind,ind+1) = 1/h^2;
        else
            % Boundary condition
            f(ind) = f(ind) - 0/h^2;
        end
        if j > 1
            A(ind,ind-N) = 1/h^2;
        else
            % Boundary condition
            f(ind) = f(ind) - 0/h^2;
        end
        if j < N
            A(ind,ind+N) = 1/h^2;
        else
            % Boundary condition
            f(ind) = f(ind) - 0/h^2;
        end
        
        % Right-hand-side
        x = x_int(i);
        y = y_int(j);
        f(ind) = f(ind) - exp(-(x-0.25)^2 - (y-0.25)^2);
    end
end
A
A\f




