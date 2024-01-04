classdef nap % Functions file
    methods(Static) % Do nap.function_name to use
        
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

          function [R_t] = cholesky(A)
            n = size(A);
            for j=1:n-1
                A(j,j) = sqrt(A(j,j));
                for i = j+1:n
                    A(i,j) = A(i,j) / A(j,j);
                end
                for k = j+1:n
                    for i = j+1:n
                        A(i,k) = A(i,k) - A(i,j)*A(k,j);
                    end
                end
            end
            A(n,n) = sqrt(A(n,n));
            R_t = tril(A);
          end
         
          function [Q, R] = gram_schmidt(A)
            [m,n] = size(A);
%             Q = zeros(m,n);
%             R = zeros(n,n);
            for j=1:n
                x=A(:,j);
                for i=1:j-1
                    R(i,j) = Q(:,i)'*A(:,j);
                    x=x-R(i,j)*Q(:,i);
                end
                R(j,j)=norm(x);
                Q(:,j) = x/R(j,j);
            end
          end

          function [Q, R] = modified_gram_schmidt(A)
            [m,n] = size(A);
%             Q = zeros(m,n);
%             Q(1:m,1) = A(1:m,1);
%             R=zeros(n);
%             R(1,1) =1;
            
            for k=1:n
                R(k,k) = norm(A(1:m,k));
                Q(1:m,k) = A(1:m,k)/R(k,k);
                for j=k+1:n
                    R(k,j) = Q(1:m,k)' * A(1:m,j);
                    A(1:m,j) = A(1:m,j) - Q(1:m,k)*R(k,j);
                end
            end
          end

          function [Q, R] = house_holder(A)
            [m, n] = size(A);
            R = A; 
            Q = eye(m,n);
            V = zeros(m,n);
            
            for k = 1:n
                x = R(k:m,k);
                
                if x(1)==0
                    s=1;
                else 
                    s=sign(x(1));
                end
                
                x(1)=s*norm(x,2)+x(1);
                
                v=x;
                v = v./norm(v,2);
                V(k:m,k)=v;
                
                R(k:m, k:n) = R(k:m, k:n) - 2*v*(v'*R(k:m, k:n));
            end
            
            for i=1:n
                for k=n:-1:1
                    Q(k:m,i) = Q(k:m,i) - 2*V(k:m,k)*(V(k:m,k)'*Q(k:m,i));
                end
            end
            R=R(1:n,:);
          end

    end
end