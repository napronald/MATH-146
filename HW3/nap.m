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
    end
end