clear all; clc; close all;

A = [1 2; 4 5];


[U,sig,V] = svd(A)


theta = linspace(0, 2*pi, 100);
x = cos(theta);
y = sin(theta);
S = [x; y];
plot(S(1,:), S(2,:), 'b');
AS = A*S;
hold on;
plot(AS(1,:), AS(2,:), 'r');
