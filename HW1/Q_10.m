clear all
clc
close all

x = [-20:1:20];

exact = exp(x);

for i = 1:length(x)
    approx(i) = myexp(x(i));
end

abs_error1 = abs(exact-approx);
rel_error1 = abs((exact-approx) ./ (exact));

subplot(2,1,1)
semilogy(x,rel_error1,'r-',x,abs_error1,'b-')
axis([-20 20 10^-17 1])
xlabel('x')
ylabel('Error')
titletxt = sprintf('Absolute and Relative Error for original version');
title(titletxt)
legend('Relative Error','Absolute Error')

for i = 1:length(x)
    approx(i) = mod_exp(x(i));
end

abs_error2 = abs(exact-approx);
rel_error2 = abs((exact-approx) ./ (exact));

subplot(2,1,2)
semilogy(x,rel_error2,'r-',x,abs_error2,'b-')
axis([-20 20 10^-27 1])
xlabel('x')
ylabel('Error')
titletxt = sprintf('Absolute and Relative Error for modified version');
title(titletxt)
legend('Relative Error','Absolute Error')


function [approx,Nterms]=mod_exp(x);
    if x > 0
        oldsum = 0;
        newsum = 1;
        term = 1;
        n = 0;
        while newsum~=oldsum
            n = n+1;
            term = term*x/n;
            oldsum = newsum;
            newsum = newsum + term;
        end
        Nterms = n + 1;
        approx = newsum;
    else
        x = abs(x);
        oldsum = 0;
        newsum = 1;
        term = 1;
        n = 0;
        while newsum~=oldsum
            n = n+1;
            term = term*x/n;
            oldsum = newsum;
            newsum = newsum + term;
        end
        Nterms = n + 1;
        approx = newsum;
        approx = 1/approx;
    end
end


function [approx,Nterms]=myexp(x);
    oldsum = 0;
    newsum = 1;
    term = 1;
    n = 0;
    while newsum~=oldsum
        n = n+1;
        term = term*x/n;
        oldsum = newsum;
        newsum = newsum + term;
    end
    Nterms = n + 1;
    approx = newsum;
end