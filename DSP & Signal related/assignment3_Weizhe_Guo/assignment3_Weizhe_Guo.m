% Name: Weizhe Guo
% Date: 02/06/2018
% Assignment 3

%% 1. Lunar Eclipse
% a)
A = ones(100);
% b)
B = zeros(100);
% c)
i = 1:100;
j = 1:100;
[I,J] = meshgrid(i,j);
M1 = sqrt((I-50).^2+(J-50).^2);
for m = 1:100
    for n =1:100
        if M1(m,n) < 20
            A(m,n) = 0;
        end
    end
end
% d)
i = 1:100;
j = 1:100;
[I,J] = meshgrid(i,j);
M2 = sqrt((I-40).^2+(J-40).^2);
for m = 1:100
    for n =1:100
        if M2(m,n) < 20
            B(m,n) = 1;
        end
    end
end
% e)
figure
imshow(A);

figure
imshow(B);

figure
imshow(A&B);

figure
imshow(A|B);

figure
imshow(~(A&B));

figure
imshow(~(A|B));

%% 2. Fun with find
Z = sin(linspace(0,5,100))+1;
[val, ind] = findClosest(Z,3/2);

%% 3. Sincing Ship
% a)
x = linspace(-2*pi, 2*pi, 10001);
y = sinc(x);
plot(x,y);
hold on;
% b)
findzero = @(x) find(x.*(circshift(x,1))<0);
% c)
xcord = x(findzero(y));
ycord = y(findzero(y));
plot (xcord,ycord,'ko');
% d)
dy = diff(y)./(4*pi/10000)
xder = x(findzero(dy));
yder = y(findzero(dy));
plot(xder,yder,'r*');

    