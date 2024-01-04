clear all; clc; close all

% u = @(x,y) -exp((y-0.25)^2 - (x-0.25)^2);

u = @(x,y) -2*(x-y)*exp((y-0.25)^2 - (x-0.25)^2);

N = 2^2;
a=0;
b=1;

h = 1 / (N + 1);
% x = a:h:b;
% y = a:h:b;

xg=h*(1:N); 
yg=h*(1:N); 

[xg,yg]=ndgrid(xg,yg);

bmat = u(xg,yg);

b = reshape(bmat,N*N,1)




% A = zeros(N*N);
% f = zeros(N*N,1);
% 
% I = speye(N,N);
% E = sparse(2:N,1:N-1,1,N,N);
% D = E+E'-2*I;
% A = kron(D,I)+kron(I,D);
% A = full(A)

% f = h^2* exp(-(x-0.25).^2 - (y-0.25).^2)'


% A = 4*eye(N*N) - diag(ones(N*N-1,1),1) - diag(ones(N*N-1,1),-1) - diag(ones(N*N-3,1),-3) - diag(ones(N*N-3,1),3)

% for i=1:N*N
%     for j=1:N*N
%         idx = (j-1)*N + i;
%         A(idx,idx) = -4;
% %         A(i,j) = 1;
% %         A(i,j) = 1;
% 
%         if i < N
%             A(idx,idx) = 1;
%         end
%         if i > 1
%             A(idx,idx) = 1;
%         end
        
        
%         else
%             % Boundary condition
%             f(ind) = f(ind) - 0/h^2;
%         end
%         A(i,i+2) = 1;
%         A(i+2,i) = 1;
%     end
% end
% A

%%
% Loop over interior grid points
% for i=1:N
%     for j=1:N
%         % Index of current grid point
%         ind = (j-1)*N + i;
%         
%         % Finite difference approximation of Laplacian
%         A(ind,ind) = -4/h^2;
%         if i > 1
%             A(ind,ind-1) = 1/h^2;
%         else
%             % Boundary condition
%             f(ind) = f(ind) - 0/h^2;
%         end
%         if i < N
%             A(ind,ind+1) = 1/h^2;
%         else
%             % Boundary condition
%             f(ind) = f(ind) - 0/h^2;
%         end
%         if j > 1
%             A(ind,ind-N) = 1/h^2;
%         else
%             % Boundary condition
%             f(ind) = f(ind) - 0/h^2;
%         end
%         if j < N
%             A(ind,ind+N) = 1/h^2;
%         else
%             % Boundary condition
%             f(ind) = f(ind) - 0/h^2;
%         end
%         
%         % Right-hand-side
%         x = x_int(i);
%         y = y_int(j);
%         f(ind) = f(ind) - exp(-(x-0.25)^2 - (y-0.25)^2);
%     end
% end
% A
% A\f





% for i=1:N
%     for j=1:N
%         ux(i,j) = (u(i+1,j) - 2*(u(i,j)) + u(i-1,j)) / h^2;
%         uy(i,j) = (u(i,j+1) - 2*(u(i,j)) + u(i,j-1)) / h^2;
%         ux + uy
%     end
% end

% f(idx) = h^2*(-exp((xg(i,j)-0.25)^2 + (yg(i,j)-0.25)^2));

%%
% N = 2^2;
% h = 1/(N+1);
% xg = h*(1:N)+h; 
% yg = h*(1:N)+h; 
% [xg, yg] = ndgrid(xg, yg);
% f = zeros(N^2, 1);
% for i = 1:N
%   for j = 1:N
%     idx = (j-1)*N + i;
%     f(idx) = h^2*(-exp((xg(i,j)-0.25)^2 + (yg(i,j)-0.25)^2));
%     if i > 1
%       f(idx) = f(idx) + 0.25;
%     end
%     if i < N
%       f(idx) = f(idx) + 0.25;
%     end
%     if j > 1
%       f(idx) = f(idx) + 0.25;
%     end
%     if j < N
%       f(idx) = f(idx) + 0.25;
%     end
%   end
% end
% f