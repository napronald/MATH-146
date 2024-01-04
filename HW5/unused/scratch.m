clear all; clc

% A = 10*randn(500,200);
% n = 10;
% A = rand(n-1,n);
% A(end+1,:) = sum(A);
% A;

% [V,R] = householder_qr(A)

% V*R

% [Q, R] = modified_gram_schmidt(A);

% Q
% R

% verify = Q*R;
% rel_error = norm(verify-A) / norm(A)
% 
% function [Q, R] = modified_gram_schmidt(A)
% [m,n] = size(A);
% 
% for i=1:n
%     R(i,i) = norm(A(1:m,i));
%     Q(1:m,i) = A(1:m,i)/R(i,i);
%     for j=i+1:n
%         R(i,j) = Q(1:m,i)' * A(1:m,j);
%         A(1:m,j) = A(1:m,j) - Q(1:m,i)*R(i,j);
%     end
% end
% end




% function [A, p] = house(A) 
% [m,n] = size(A);
% 
% R = A;
% for k = 1:m-1
%     x(k:m,1) = R(k:m, k);
%     g = norm(x);
%     v = x;
%     v(k) = x(k) + g;
%     s = norm(v);
%     if s ~= 0
%         w = v/s;
%         R = R-w*u';
%         Q = Q-2*Q*W*W';
%     end
% end
% end


% function [Q, R] = modified_gram_schmidt(A)
% [m,n] = size(A);
% Q = zeros(m,n);
% R=zeros(n);
% 
% % Q(1:m,1) = A(1:m,1);
% % R(1,1) = 1;
% 
% V=A;
% 
% for i=1:n
%     R(i,i) = norm(V(1:m,i));
%     Q(1:m,i) = V(1:m,i)/R(i,i);
%     for j=i+1:n
%         R(i,j) = Q(1:m,i)' * V(1:m,j);
%         v(1:m,j) = V(1:m,j) - Q(1:m,i)*R(i,j);
%     end
% end
% end


% function [V,R] = householder_qr(A)
% [m,n] = size(A);
% Q = eye(m);
% for k=1:min(m-1,n)
% ak = A(k:m,k); % vector to be zeroed out
% vk = ak; vk(1) = vk(1) + sign(ak(1))*norm(ak); % 1. construct vk that defines the reflector
% vk = vk/norm(vk); % 2. normalize vk
% A(k:m,k:n) = A(k:m,k:n) - 2*vk*(vk*A(k:m,k:n)); % 3. update A
% V = vk; % store vk in a cell array
% end
% R = A;
% end



% function qr_reorth_rank(matrix, reorth_threshold=10, rank_deficiency_factor=100, update_r=false)
%     num_vectors = size(matrix)[2]
%     orth_matrix = copy(matrix)
%     r_matrix = zeros(num_vectors, num_vectors)
%     for vec_idx = 1:num_vectors
%         norm_orth_column = 0
%         perform_reorthogonalisation = true
%         while (perform_reorthogonalisation)
%             norm_current_column = norm(orth_matrix[:, vec_idx])
%             for span_base_idx = 1:(vec_idx - 1)
%                 projection_length = dot(orth_matrix[:, span_base_idx], orth_matrix[:, vec_idx])
%                 if (norm_orth_column == 0)
%                     r_matrix[span_base_idx, vec_idx] = projection_length
%                 elseif (update_r)
%                     r_matrix[span_base_idx, vec_idx] += projection_length
%                 end
%                 orth_matrix[:, vec_idx] -= projection_length * orth_matrix[:, span_base_idx]
%             end
%             norm_orth_column = norm(orth_matrix[:, vec_idx])
%             if ((norm_orth_column < norm_current_column / reorth_threshold) &&
%                 (norm_orth_column > rank_deficiency_factor * eps() * norm_current_column))
%                 perform_reorthogonalisation = true
%             else
%                 perform_reorthogonalisation = false
%                 if (norm_orth_column <= rank_deficiency_factor * eps() * norm_current_column)
%                     orth_matrix[:, vec_idx] .= 0
%                 end
%             end
%         end
%         if (norm(orth_matrix[:, vec_idx]) > eps())
%             orth_matrix[:, vec_idx] = orth_matrix[:, vec_idx] / norm(orth_matrix[:, vec_idx])
%         else
%             orth_matrix[:, vec_idx] .= 0
%         end
%         r_matrix[vec_idx, vec_idx] = dot(orth_matrix[:, vec_idx], matrix[:, vec_idx])
%     end
%     return (orth_matrix, r_matrix)
% end


% function [Q, R] = modified_gram_schmidt(A)
% [m,n] = size(A);
% if m < n
%     error('m < n')
% end
% 
% Q = zeros(m,n);
% Q(1:m,1) = A(1:m,1);
% R=zeros(n);
% R(1,1) =1;
% 
% for k=1:n
%     R(k,k) = norm(A(1:m,k));
% %     if R(k,k) == 0
% %         Q(1:m:k) = 0;
% %     else
%         Q(1:m,k) = A(1:m,k)/R(k,k);
% %     end
%     for j=k+1:n
% %         if R(k,k) == 0
% %             R(k,j) = A(1:m,j);
% %             A(1:m,j) = A(1:m,j) - R(k,j);
% %         else
%             R(k,j) = Q(1:m,k)' * A(1:m,j);
%             A(1:m,j) = A(1:m,j) - Q(1:m,k)*R(k,j);
% %         end
%     end
% end
% end

% function [Q, R] = modified_gram_schmidt(A)
% [m,n] = size(A);
% 
% % Q = zeros(m,n);
% % Q(1:m,1) = A(1:m,1);
% % R=zeros(n);
% % R(1,1) =1;
% 
% for i=1:n
%     R(i,i) = norm(A(1:m,i));
%     Q(1:m,i) = A(1:m,i)/R(i,i);
%     for j=i+1:n
%         R(i,j) = Q(1:m,i)' * A(1:m,j);
%         A(1:m,j) = A(1:m,j) - Q(1:m,i)*R(i,j);
%     end
% end
% end
%%
clear all; clc;


A = [1,2;1,4]

[W,R] = house(A)
W*R

function [W,R] = house(A)
[m,n] = size(A);
W = zeros(m,n);
y = eye(m,m);

for k=1:n
    x = A(k:m,k);
    e1 = [1; zeros(m-k,1)];
    if x(1) >= 0
        v = x + norm(x)*e1;
    else
        v = x-sign(x(1))*norm(x)*e1;
    end
    v = v/norm(v);
    A(k:m,k:n) = A(k:m,k:n)-2*v*(v'*A(k:m,k:n));
    W(k:m,k) = v;
end
R = A(1:n,1:n);
end


%% 
clear all; clc;

% n =4;
% A = rand(n-1,3);
% A(end+1,:) = sum(A);
% A

% A = 4x3 to 4x2 solution 2x1

A = 10*randn(4,3);
A(:,2) = 2*A(:,1) ;
A

b = randn(4,1)

[Q,R] = gram_schmidt(A)

A(:,2) = []

A \ b

y = Q'*b;

coe = nap.back_sub(R,y)
% AA = Q*R

% M = AA'*AA;
% y = AA'*b;

% R_t = nap.cholesky(M);
% R = transpose(R_t);

% [z] = nap.forward_sub(R_t,y)
% [coefs] = nap.back_sub(R,z)
% 
% AA(:,n) = mat
% c = Q'*b
% [x] = nap.back_sub(R,c)

% verify = A \ b
% rel_error = norm(x-verify) / norm(verify)

function [Q, R] = gram_schmidt(A)
tol=1e-10;
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);
n_discard = 0;
for j=1:n 
    x=A(:,j);
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        x=x-R(i,j)*Q(:,i);
    end
%     norm(x);
    R(j,j)=norm(x);
    if norm(x) < tol
        n_discard = n_discard + 1;
        discard(n_discard) = j;
%         R(:,j) = [];
%         R(j,:) = [];
%         Q(:,j) = [];
    else
        Q(:,j) = x/R(j,j);
    end
    
end
for n=1:length(discard)
    R(:,discard) = [];
    R(discard,:) = [];
    Q(:,discard) = [];
end
end


%% working code for the bonus


clear all; clc;

% A = 4x3 to 4x2 solution 2x1

A = 10*randn(40,30);
A(:,2) = 2*A(:,1) ;
% A(:,5) = 2*A(:,4);
A
b = randn(40,1)

[Q,R] = gram_schmidt(A)
% size(R)

y = Q'*b;

coe = nap.back_sub(R,y)

% size(A)
A(:,2) = [];
% A(:,5) = [];
% size(A)

coe_matlab = A \ b


rel_diff = norm(coe-coe_matlab,2) / norm(coe_matlab)


function [Q, R] = gram_schmidt(A)
tol=1e-12;
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);
n_discard = 0;
for j=1:n 
    x=A(:,j);
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        x=x-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(x);
    if norm(x) < tol
        n_discard = n_discard + 1;
        discard(n_discard) = j;
    else
        Q(:,j) = x/R(j,j);
    end
    
end
% for n=1:length(discard)
R(:,discard) = [];
R(discard,:) = [];
Q(:,discard) = [];
% end
end











