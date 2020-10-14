% Name: Weizhe Guo
% Date: 02/13/2018
% Assignment 4

% Acknowledgement: Thanks Yilun Eric Jiang for telling me the
% convinient method to check ONS. Thanks Zheng Alex Liu for explaining the
% last problem and ndgrid

% d)
testvec = complex(rand(4,1),rand(4,1));
orthset1 = gramschmidt(complex(rand(4,3),rand(4,3))); 
orthset2 = gramschmidt(complex(rand(4,4),rand(4,4)));
onstest1 = isOrthnormal(orthset1);
onstest2 = isOrthnormal(orthset2);
orthproj1 = orthoProj(testvec,orthset1);
orthproj2 = orthoProj(testvec,orthset2);
error1 = testvec - orthproj1
error2 = testvec - orthproj2
% error2 is much smaller because the dimension of the vector matches the
% dimension of the O.N.S. So we actually generate an O.N.B.

[x,y] = ndgrid(0:0.01:2*pi,0:pi/2:2*pi);
z1 = sin(x);
z2 = 1/sqrt(2*pi)*exp(-(x-y).^2);
plot(x,z1)
hold on
plot(x,z2)
title('Plot of sine and Gaussians');
xlabel('x');
ylabel('y');
axis([0 2*pi -1 1]);
legend('sine','Gaussian1','Gaussian2','Gaussian3','Gaussian4','Gaussian5')

orthGaussian = gramschmidt(z2);
sinproj = orthoProj(z1,orthGaussian);
figure
subplot(2,1,1)
plot(x,sinproj)
hold on
plot(x,z1)
title('Sine and its Projection on ONS');
xlabel('x');
ylabel('y');
legend('Projection','Sin');

subplot(2,1,2)
plot(x,orthGaussian)
hold on
title('O.N.B. functions');
xlabel('x');
ylabel('y');
legend('ON Gaussian1','ON Gaussian2','ON Gaussian3','ON Gaussian4','ON Gaussian5');