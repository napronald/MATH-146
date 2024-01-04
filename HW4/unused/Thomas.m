a = [0 1 1 -1];
b = [2 -1 1 1];
c = [6 2 2 0];
d = [14; -3; 9; 5];

x = thomas(a,b,c,d)

function [x] = thomas(a, b, c, d)
%function to perform Thomas Algorithm to solve tridiagonal system with
%vectors a,b,c on the subdiagonal, diagonal, and superdiagonal, resp.
%a,b,c nx1 first term of a not used; last term of c not used
%d, rhs
%output: x to solve Ax=b

n=length(d);
if length(a)~=n || length(b)~=n || length(c)~=n
    error('Incompatible vector lengths')
end

cp=zeros(n,1);  %cprime (from class)
dp=zeros(n,1);  %dprime (from class)
x=zeros(n,1);

dpdenom=b(1);   %denominator of dprime (from class)
dp(1)=d(1)/dpdenom;

for i=2:n
    cp(i-1)=c(i-1)/dpdenom;     %fill in cprime
    dpdenom=b(i)-a(i)*cp(i-1);  %make new dprime denominator
    dp(i)=(d(i)-a(i)*dp(i-1))/dpdenom;   %dprime
end

x(n)=dp(n);
for i=n-1:-1:1
    x(i)=dp(i)-cp(i)*x(i+1);
end

end




