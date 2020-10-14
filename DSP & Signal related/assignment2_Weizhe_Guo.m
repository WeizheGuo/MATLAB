% Name: Weizhe Guo
% Date: 02/04/2018
% Assignment 2

%% 1. Vector
% a)
x = linspace(0,1,100);
y = exp(x);
% b)
s1 = trapz(y)*1/99;
s2 = sum(y)*1/99; 
% c)
s3 = cumtrapz(x)*1/99;
s4 = cumsum(x)*1/99;
% d)
dx = diff(x)/(1/99);
dx2 = diff(dx)/(1/98);
% The size of each derivative is one less than the previous vector
% e)
vec1 = circshift([0 1 2 3 4 5],3);

%% 2. Array Foray
% a)
A = reshape(1:1:100,[10 10]);
% b)
B = magic(10);
C = diag(diag(B));
% c)
B(:,2) = flipud(B(:,2));
% d)
A = fliplr(A);
% e)
csum = sum(A*B,1);
% f)
cmean = sum(A.*B,2)/10;
% g)
A(:,10)=[] 

%%3. Loops and time
% a)
tic
for i = 1:300
    for j = 1:500
        a1 = (i*i+j*j)/(i+j+3);
    end
end
toc
t1 = toc;
% b)
a2 = zeros(300,500);
tic
for i = 1:300
    for j = 1:500
        a2 = (i*i+j*j)/(i+j+3);
    end
end
toc
t2 = toc;
% c)
tic
m = 1:500;
n = 1:300;
[M,N]= meshgrid(m,n);
a3 = (M.^2+N.^2)./(M+N+3);
toc
t3 = toc;

table = table({'T1';'T2';'T3'}, [t1; t2; t3]);
