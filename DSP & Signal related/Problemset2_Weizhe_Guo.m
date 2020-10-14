% Weizhe Guo
% ECE-211
% 02/07/2018
% Problem Set 1

h = [2 -1 2 3 -1];
x = [-2 4 1 1];
y = conv(h,x);
n = [-3:1:4];
stem(n,y)