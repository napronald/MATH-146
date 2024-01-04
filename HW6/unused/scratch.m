clear all
clc
close all
% Grid size
N = 3;

% Grid spacing
h = 1/(N+1);

% Define the grid
xg = h*(1:N);
yg = h*(1:N);
[xg,yg] = ndgrid(xg,yg);

% Interior grid points
x_int = xg + h;
y_int = yg + h;

% Initialize Laplacian matrix and right-hand-side vector
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
% u = @(x,y) exp(x*y)+x^2;
% x = 1;
% y = 2;
% h = 2.^(-2:-1:-12);
% delta_u = zeros(size(h));
% for i = 1:length(h)
%     delta_u(i) = laplacian_2d(u, x, y, h(i));
% end
% 
% 
% exact_delta_u = 2*exp(2);
% rel_err = abs(delta_u - exact_delta_u) ./ abs(exact_delta_u);
% 
% 
% 
% figure;
% loglog(h, h.^2, '-o');
% xlabel('h');
% ylabel('Relative error');
% legend('Computed relative error', 'h^2');
% 
% 
% % function delta_u = laplacian_2d(u, x, y, h)
% %     delta_u = (u(x-h,y) + u(x+h,y) + u(x,y-h) + u(x,y+h) - 4*u(x,y)) / h^2;
% % end
% 
% %%
% % Define the PDE
% u = @(x) (-1/pi)*x.*cos(pi*x) + 1/pi + 1;
% 
% % Set up the grid
% N_values = 2.^(3:12);
% err = zeros(size(N_values));
% for i = 1:length(N_values)
%     N = N_values(i);
%     x = linspace(0,1,N+2)';
%     h = 1/(N+1);
% 
%     % Construct the matrix and right-hand side
%     A = (2/h^2)*eye(N) - (1/h^2)*diag(ones(N-1,1),1) - (1/h^2)*diag(ones(N-1,1),-1);
%     f = h^2*sin(pi*x(2:end-1));
% 
%     % Solve the system
%     u_approx = [0; A\f; 1];
% 
%     % Compute the error
%     err(i) = max(abs(u(x) - u_approx));
% end
% 
% % Plot the error vs N on a log-log scale
% figure
% loglog(N_values, err, 'o-', 'LineWidth', 2)
% xlabel('N')
% ylabel('Relative error in max norm')
% 
% % Plot the true and approximate solutions for the finest grid
% x_fine = linspace(0,1,1000)';
% u_fine = u(x_fine);
% N_fine = N_values(end);
% x = linspace(0,1,N_fine+2)';
% h = 1/(N_fine+1);
% A = (2/h^2)*eye(N_fine) - (1/h^2)*diag(ones(N_fine-1,1),1);

%%
clear all

global num_cell N num_elements

%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503*10^-23;          %Boltzmann const., J/k
T = 296.;                      %temperature
epsilon_0 =  8.85418782*10^-12; %F/m
Vt = (kb*T)/q;

%% System Setup
tic

num_cell = 1001;
N = num_cell -1;   %number of INTERIOR mesh points
num_elements = N^2;  %NOTE: this will specify number of elements in the solution vector V which = (num_cell +1 -2)^2 b/c we are in 2D
%(num_cell +1) = # of mesh pts, b/c matlab starts indexing from 1, then -2
%the endpts

Va = 1.; %applied voltage

%Matrix of the system's net charge
%for now make the net charge = 0  --> this will give a linear 2D potential
%plot
netcharge = zeros(num_cell+2,num_cell+2);
%netcharge = ones(num_elements,num_elements);  %this will make a slightly curved 2D potential plot

%% Define dielectric constant matrix
%NOTE: epsilons will be defined at 1/2 integer points, so epsilons inside
%the cells, not at cell boundaries
%will use indexing: i + 1/2 is defined as i+1 for the index
epsilon = zeros(num_cell+2, num_cell +2);
%later will fill with real epsilons corresponding to the different layers

%NOTE: I unfortunately can't define epsilon(0,..) in matlab, so the endpts
%are at 1...
for i = 1:num_cell+2
    epsilon(i,:) = 1;
end

%% Define boundary conditions and initial conditions
V_bottomBC = 0;
V_topBC = Va/Vt;

% Initial conditions
diff = (V_topBC - V_bottomBC)/num_cell;
index = 0;
for j = 1:N  %corresponds to z coord
    index = index +1;
    V(index) = diff*j;
    for i = 2:N  %elements along the x direction assumed to have same V
        index = index +1;
        
        V(index) = V(index-1);
    end
end

%Side BCs will be filled in from V, since are insulating BC's
%i's in these BC's correspond to the x-value (z values are along a line,
%top and bottom)
%THESE NEED TO BE UPDATED AT EVERY ITERATION OF POISSON SOLVE
for i = 1:N
    V_leftBC(i) = V((i-1)*N + 1);  %just corresponds to values of V in the 1st subblock
    V_rightBC(i) = V(i*N);
end


%% Set up matrix equation and solve
AV = SetAV_2D(epsilon);
%spy(AV);  %allows to see matrix structure, very useful!

%set up rhs of Poisson equation. Note for epsilons, are assuming that
%epsilons at the boundaries are the same as espilon 1 cell into interior of
%device
index = 0;
for j = 1:N
    if(j ==1)  %different for 1st subblock
        for i = 1:N
            index = index +1;
            if (i==1)  %1st element has 2 BC's
                bV(index,1) = netcharge(i+1,j+1) + epsilon(i+1,j+1)*(V_leftBC(1) + V_bottomBC);  %+1 b/c netcharge and epsilon include endpoints but i,j index only the interior
            elseif (i==N)
                bV(index,1) = netcharge(i+1,j+1) + epsilon(+1,j+1)*(V_rightBC(1) + V_bottomBC);
            else
                bV(index,1) = netcharge(i+1,j+1) + epsilon(i+1,j+1)*V_bottomBC;
            end
        end
    elseif(j == N)  %different for last subblock
        for i = 1:N
            index = index +1;
            if (i==1)  %1st element has 2 BC's
                bV(index,1) = netcharge(i+1,j+1) + epsilon(i+1,j+1)*(V_leftBC(N) + V_topBC);
            elseif (i==N)
                bV(index,1) = netcharge(i+1,j+1) + epsilon(i+1,j+1)*(V_rightBC(N) + V_topBC);
            else
                bV(index,1) = netcharge(i+1,j+1) + epsilon(i+1,j+1)*V_topBC;
            end
        end
    else %interior subblocks
        for i = 1:N
            index = index +1;
            if(i==1)
                bV(index,1) = netcharge(i+1,j+1) + epsilon(i+1,j+1)*V_leftBC(j);
            elseif(i==N)
                bV(index,1) = netcharge(i+1,j+1) + epsilon(i+1,j+1)*V_rightBC(j);
            else
                bV(index,1) = netcharge(i+1,j+1);
            end
        end
    end
end

%solve for V
%spparms('spumoni',2)
V = AV\bV;

toc

%plot V
surf(1:N,1:N,reshape(V,N,N))  %need to reshape the vector into a 2D matrix, to make a 2D plot


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