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