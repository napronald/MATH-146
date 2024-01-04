clear all
close all
clc

% t =[0 1 2];
% b = [0.1 0.9 2.0];

% t = [-2 -1 0 1 2 3];
% b = [9;5;3;4;8;12];
% n = 3;
% coefs = lsfit(t,b,n)

tt = [-2 -1 0 1 2 3];
bb = [9;5;3;4;8;12];
terms = 5;


% m = 21;
% tt = 0:1/(m-1):1;
% bb = cos(2*pi*tt);

for n=1:terms
    coefs{n} = lsfit(tt,bb,n);
end

t = min(tt):0.55:max(tt);
length(t)
z = ones(terms,length(t));

for n = 1:terms
    z(n,:) = z(n,:).*t + coefs{n}(n);
    for j=n-1:-1:1
        z(n,:) = z(n,:).*t + coefs{n}(j);
    end
end
z
plot(tt,bb,'ro')
hold on
plot(t,z)
xlabel('t')
ylabel('p_{n-1}')
legend('real','1','2','3','4','5')

function coefs = lsfit(t,b,n)
t = t(:);
b = b(:);
m = length(t);

A = ones(m,n);

for j=1:n-1
    A(:,j+1) = A(:,j).*t;
end

B = A'*A;
y = A'*b;

coefs = B \ y;
end