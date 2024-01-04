clear all; clc; close all;

load mandrill.mat
X1 = X;
[m1,n1] = size(X1);

% load durer.mat
% X2 = X;
% [m2,n2] = size(X2);

% Compute SVD and truncate for various ranks
r_values = [2 4 8 16 32 64 128 256];

% Mandrill image
[U1,S1,V1] = svd(X1);
figure;
for i = 1:length(r_values)
    r = r_values(i);
    X1_truncated = U1(:,1:r)*S1(1:r,1:r)*V1(:,1:r)';
    subplot(2,4,i);
    imagesc(X1_truncated);
    colormap(gray);
    title(sprintf('Rank = %d',r));
    axis off;
end

% 
% % Durer image
% [U2,S2,V2] = svd(X2);
% figure;
% for i = 1:length(r_values)
%     r = r_values(i);
%     X2_truncated = U2(:,1:r)*S2(1:r,1:r)*V2(:,1:r)';
%     subplot(2,4,i);
%     imagesc(X2_truncated);
%     colormap(gray);
%     title(sprintf('Rank = %d',r));
%     axis off;
% end


