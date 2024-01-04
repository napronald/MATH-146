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
            V=A;
            for i=1:n
                R(i,i) = norm(V(1:m,i));  
                Q(1:m,i) = V(1:m,i)/R(i,i);
                for j=i+1:n
                    R(i,j) = Q(1:m,i)' * V(1:m,j); 
                    V(1:m,j) = V(1:m,j) - Q(1:m,i)*R(i,j);
                end
            end
          end

          function [Q,R] = house_holder(A)
            [m, n] = size(A);
            R = A; 
            Q = eye(m,n);
            for k = 1:n
                x = R(k:m,k);
                if x(1) == 0 
                    x(1)= norm(x,2) + x(1);
                else 
                    x(1)= sign(x(1))*norm(x,2) + x(1);
                end
                x = x/norm(x,2); 
                V(k:m,k)= x/norm(x,2);
                R(k:m, k:n) = R(k:m, k:n) - 2*x*(x'*R(k:m, k:n)); 
            end
            for j=1:n
                for i=n:-1:1 
                    Q(i:m,j) = Q(i:m,j) - 2*V(i:m,i)*(V(i:m,i)'*Q(i:m,j));
                end
            end
            R=R(1:n,:);
          end

    end
end