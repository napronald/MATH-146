classdef nap % Functions file
    methods(Static) % Do nap.function_name to use
          function [L, U] = LU_decomp(A)
            n = size(A, 1); 
            L = eye(n); % Initialize L,U
            U = A;
            for j = 1:n-1 % loop over columns,rows
                for i=j+1:n
                    % Assign values to L,U components
                    L(i,j) = U(i,j) / U(j,j); 
                    U(i,j) = 0;
                    U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
                end
            end
            D = diag(U); 
            search = find(D==0); % Check for 0 in diagonal
            if ~isempty(search)
                error('Zero found in diagonal')
            end
          end
          
          function [x] = forward_sub(L,B)
            n = length(B);
            x = zeros(length(B),1); % Initialize solution vector
            
            x(1) = B(1)/L(1,1); % Find first entry
            for i = 2:n 
                sum = 0;
                for j = 1:i-1 
                    sum = sum + L(i,j) * x(j); % Compute sum
                end
                x(i) = (B(i) - sum) / L(i,i); % Compute x entries
            end
          end

          function [x] = back_sub(U,B)
            n = length(B);
            x = zeros(length(B),1); % Initialize solution vector
            
            x(n) = B(n)/U(n,n); % find last entry
            for i = n-1:-1:1 
                sum = 0;
                for j = n:-1:1 
                    sum = sum + U(i,j) * x(j); % Compute sum 
                end
                x(i) = (B(i) - sum) / U(i,i); % Compute x entries
            end
          end

          function [P, L, U] = LUP(A)
            n=length(A); % Initialize P,L,U
            L=eye(n); 
            P=eye(n);
            U=A;
            for j=1:n-1 % loop over columns
                [max_num, index] = max(abs(U(j:n,j))); % find index 
                index = index+j-1; % adjust indices based on size of column
                if index ~= j % Apply row swaps
                    U([index, j], :) = U([j, index], :);
                    P([index, j], :) = P([j, index], :);
                    L([index, j], 1:j-1) = L([j, index], 1:j-1);
                end
                for i=j+1:n % Assign values to L,U components
                    L(i,j) = U(i,j) / U(j,j);
                    U(i,j) = 0;
                    U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
                end
            end
          end

          function [x] = Ax_B(A,B)
            % Using functions from nap.m file
            [P, L, U] = nap.LUP(A);
            
            B = P*B; % formulate PB as B
            
            % Apply Ly=B --> Ux=y
            y = nap.forward_sub(L,B);
            x = nap.back_sub(U,y);
          end
          
          function [A,b] = rand_sdd_matrix(n)
            rng('shuffle');
            A=rand(n);
            B=A+A';% or A*A'
            C=diag(diag(B)); % extraction of the diagonal part of B
            D=B-C;
            m=max(sum(D)); % m is the biggest sum in columns (or in rows) of D 
%             E=D+m*eye(n);% the zero diagonal of D is filled with values m
            E=D+(m+0.001)*eye(n);
            A = E;
            b = rand(n,1);
          end

          
          
          function [L] = lap1d(n)
            %Creates the 1-D discretized laplacian operator (without the 1/h^2 factor)
            e=ones(n,1);
            L=spdiags([e -2*e e], [-1 0 1], n, n);      
          end

          function [L2] = lap2d(nx,ny)
            %Creates the 2-D discrete laplacian matrix (without the 1/h^2 factor)
            Lx=nap.lap1d(nx);
            Ly=nap.lap1d(ny);
            
            Ix=speye(nx);
            Iy=speye(ny);
            
            L2=kron(Iy,Lx)+kron(Ly,Ix);
            
          end

          function [u,n_it] = FDJacobi_2D(f, tol,N, xmin,ymin, L )
            %2D Jacobi method for solving 2D Poisson equation, Lapl u = f 
            % with homogeneous Dirichlet boundary conditions
            
            %Inputs: f, a function that gives right hand side
            %    tolerance
            %    N, number of interior gridpoints in each direction
            %    xmin, ymin for domain and L, length of domain (assumed square domain)
            
            %outputs: solution u, as NxN matrix
            %         n_int, number of iterations required
            
            u=zeros(N,N); %initialize
            
            h=L/(N+1);
            xg=h*(1:N)+xmin;
            yg=h*(1:N)+ymin;
            [xg,yg]=ndgrid(xg,yg);
            
            rhs=f(xg,yg);
            %make a larger matrix with 0s to deal with the boundary terms (boundary
            %terms will not be used from this matrix)
            rhstouse=zeros(N+2,N+2);
            rhstouse(2:N+1,2:N+1)=rhs;
            
            reldiff2norm=1; %initialize for stopping
            n_it=0;
            
            while reldiff2norm >tol
                %make a larger matrix with 0s to deal with the boundary terms:
                utouse=zeros(N+2,N+2);
                utouse(2:N+1,2:N+1)=u;
            
                v=zeros(size(utouse));
                for i=2:N+1
                    for j=2:N+1
            
                        v(i,j)=1/4*(utouse(i-1,j)+utouse(i,j-1)+utouse(i+1,j)+...
                            utouse(i,j+1))-h^2/4*rhstouse(i,j);
                    end
                end
                v=v(2:N+1,2:N+1);
                reldiff2norm=norm(u-v,2)/norm(u,2);
                n_it=n_it+1;
                u=v;
            end
          end

          function [u,n_it] = FDGaussSeidel_2D(f, tol,N, xmin,ymin, L )
            %2D Jacobi method for solving 2D Poisson equation, Lapl u = f 
            % with homogeneous Dirichlet boundary conditions
            
            %Inputs: f, a function that gives right hand side
            %    tolerance
            %    N, number of interior gridpoints in each direction
            %    xmin, ymin for domain and L, length of domain (assumed square domain)
            
            %outputs: solution u, as NxN matrix
            %         n_int, number of iterations required
            
            u=zeros(N,N); %initialize
            
            h=L/(N+1);
            xg=h*(1:N)+xmin;
            yg=h*(1:N)+ymin;
            [xg,yg]=ndgrid(xg,yg);
            
            rhs=f(xg,yg);
            %make a larger matrix with 0s to deal with the boundary terms (boundary
            %terms will not be used from this matrix)
            rhstouse=zeros(N+2,N+2);
            rhstouse(2:N+1,2:N+1)=rhs;
            
            reldiff2norm=1; %initialize for stopping
            n_it=0;
            
            while reldiff2norm >tol
                %make a larger matrix with 0s to deal with the boundary terms:
                utouse=zeros(N+2,N+2);
                utouse(2:N+1,2:N+1)=u;
            
                v=zeros(size(utouse));
                for i=2:N+1
                    for j=2:N+1
            
                        v(i,j)=1/4*(v(i-1,j)+v(i,j-1)+utouse(i+1,j)+...
                            utouse(i,j+1))-h^2/4*rhstouse(i,j);
                    end
                end
                v=v(2:N+1,2:N+1);
                reldiff2norm=norm(u-v,2)/norm(u,2);
                n_it=n_it+1;
                u=v;
            end
          end

          function [u,n_it] = FDSOR_2D(f, tol,N, xmin,ymin, L )
            %2D Jacobi method for solving 2D Poisson equation, Lapl u = f 
            % with homogeneous Dirichlet boundary conditions
            
            %Inputs: f, a function that gives right hand side
            %    tolerance
            %    N, number of interior gridpoints in each direction
            %    xmin, ymin for domain and L, length of domain (assumed square domain)
            
            %outputs: solution u, as NxN matrix
            %         n_int, number of iterations required
            
            u=zeros(N,N); %initialize
            
            h=L/(N+1);
            xg=h*(1:N)+xmin;
            yg=h*(1:N)+ymin;
            [xg,yg]=ndgrid(xg,yg);
            
            rhs=f(xg,yg);
            %make a larger matrix with 0s to deal with the boundary terms (boundary
            %terms will not be used from this matrix)
            rhstouse=zeros(N+2,N+2);
            rhstouse(2:N+1,2:N+1)=rhs;
            
            reldiff2norm=1; %initialize for stopping
            n_it=0;
            w=2-2*pi*h;
            while reldiff2norm >tol
                %make a larger matrix with 0s to deal with the boundary terms:
                utouse=zeros(N+2,N+2);
                utouse(2:N+1,2:N+1)=u;
            
                v=zeros(size(utouse));
                for i=2:N+1
                    for j=2:N+1
            
                        v(i,j)=w/4*(v(i-1,j)+v(i,j-1)+utouse(i+1,j)+...
                            utouse(i,j+1))-h^2*w/4*rhstouse(i,j)+(1-w)*utouse(i,j);
                    end
                end
                v=v(2:N+1,2:N+1);
                reldiff2norm=norm(u-v,2)/norm(u,2);
                n_it=n_it+1;
                u=v;
            end
          end

    end
end