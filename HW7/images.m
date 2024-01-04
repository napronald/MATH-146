clear all; clc; close all;

load mandrill.mat

[U,S,V] = svd(X);

for i = 1:8
    r = 2^i;
    A = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    subplot(2,4,i);
    imagesc(A);
    colormap(gray);
    title(['Rank = ', num2str(r)]);
end


%% 
clear all; clc; close all;

load durer.mat

[U,S,V] = svd(X);

for i = 1:8
    r = 2^i;
    A = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    subplot(2,4,i);
    imagesc(A);
    colormap(gray);
    title(['Rank = ', num2str(r)]);
end
