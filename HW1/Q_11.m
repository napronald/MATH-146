clear all
clc
close all

x = [-100:1:100];
exact = exp(x);

terms = 12;

for i = 1:length(x)
    approx(i) = myexp(x(i),terms);
end


abs_error = abs(exact-approx);
rel_error = abs((exact-approx)./(exact));
 

subplot(2,1,1)
semilogy(x,rel_error,'r-',x,abs_error,'b-')
axis([-100 100 10^-16 10^5])
xlabel('x')
ylabel('Error')
titletxt = sprintf('Absolute and Relative Error for n = 12');
title(titletxt)
legend('Relative Error','Absolute Error')


terms = 50;
for i = 1:length(x)
    approx(i) = myexp(x(i),terms);
end


abs_error = abs(exact-approx);
rel_error = abs((exact-approx)./(exact));


subplot(2,1,2)
semilogy(x,rel_error,'r-',x,abs_error,'b-')
axis([-100 100 10^-27 10^1])
xlabel('x')
ylabel('Error')
titletxt = sprintf('Absolute and Relative Error for n = 50');
title(titletxt)
legend('Relative Error','Absolute Error')



function [approx,Nterms]=myexp(x,terms);
    if x > 0
    newsum = 1;
    term = 1;
    n = 0;
    while n ~= terms
        n = n+1;
        term = term*x/n;
        newsum = newsum + term;
    end
    Nterms = n;
    approx = newsum;
    else
        x = abs(x);
        newsum = 1;
        term = 1;
        n = 0;
        while n ~= terms
            n = n+1;
            term = term*x/n;
            newsum = newsum + term;
        end
        Nterms = n;
        approx = newsum;
        approx = 1/approx;
    end
end