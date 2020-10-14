% Name: Weizhe Guo
% Date: 01/25/2018
% Assignment 1

%% 1. Creating Scalar Variables
% a)
a = 5.7*pi/6.9;
% b)
b = 239+exp(5)-2.5*10*10.^23;
% c)
c = log(4.23)*asin(0.7);
% d) 
z = (3+2j)*(4+5j);

%% 2. Complex Operations
realpart = real(z);
imagpart = imag(z);
magnitude = abs(z);
phase = angle(z);
conjugate = conj(z);

%% 3. Vector and Matrix Variable
% a)
aVec = [3.14 15 9 26+0.1j];
A1 = repmat (aVec,3,0);
A2 = A1;
% b)
bVec1 = [3.14; 15; 9; 26+0.1j];
bVec2 = aVec.';
% c)
cVec = -5:0.1:5;
% d)
dVec = transpose (linspace(-5,5,100));
% e)
B = [1+2j, 10.^(-5); exp(1j*2*pi), 3+4j];
% f)
% e)
speye(1000000);

%% 4. Vector and Matrix Operations
% a)
A = magic(5)/65;
% b)
B = randn([5 5]);
% c)
C = A*B';
% d)
D = A.*B;
% e)
F = 1/4*A.^3 + 1/4*A.^2 + 1/3*A + 1/6*eye(5);
% f)
G = inv (A);