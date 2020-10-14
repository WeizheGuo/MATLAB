% Weizhe Guo
% ECE 310
% Problem Set 3

%% Problem 2
% b)
w1= linspace(-pi,pi,100);
w2= linspace(-pi,pi,100);
h = 1/6*[1,4,1;4,-20,4;1,4,1];
Hw = freqz2(h,w1,w2);
Hwimagmax = max(abs(imag(Hw(:))));
Hwreal = real(Hw);
[ww1,ww2] = meshgrid(w1,w2);
myHw = 1/6*(-20 + 8*cos(ww1) + 8*cos(ww2) + 2*cos(ww1+ww2) + 2*cos(ww1-ww2));
maxDiff = max(abs(Hwreal-myHw));

% c)
figure
contour(w1,w2,Hw)
hold on
axis('equal')
xlabel('w1')
ylabel('w2')
title('Contour plot of H')

figure
surf(w1,w2,Hw)
axis('equal')
xlabel('w1')
ylabel('w2')
title('Surface plot of H')

% It appears to be a bandpass filter

% d)4
w3= linspace(0,0,100);
H4 = freqz2(h,w1,w3);
figure();
plot(w1,H4(1,:))
hold on
plot(w1,w1.^2)
xlabel('w')
ylabel('Frequency Response')
legend("H", "Laplacian")
title('Comparison between H and actual Laplacian')

% e)
[Hf, f1, f2] = freqz2(h);
[ff1,ff2] = meshgrid(f1,f2);
ApprLap = -(ff1.^2 + ff2.^2) * pi;
figure;
surf(f1,f2,ApprLap);
hold on;
surf(f1,f2,Hf);
hold on;
xlabel('f1')
ylabel('f2')
hold off
title('Superimposed graph of actual Lapacian and H')

%% Problem 3
% a) is in a seperate file

% b)
Hexp = InsertZeros(h);
[Hfexp,f1,f2] = freqz2(Hexp);
figure();
surf(f1,f2,Hfexp);
xlabel('f1')
ylabel('f2')
title('Contour plot for G')
figure();
contour(f1,f2,Hfexp)
xlabel('f1')
ylabel('f2')
title('Surface plot for G')
%It is a plausible approximation, but a lowpass filter might be needed. The
%upsampling may cause aliasing and thus we need an anti-aliasing filter,
%which is a lowpass filter.

%% Problem 4
% a)
hy = 1/8*[-1,-2,-1;0,0,0;1,2,1];
hx = hy';
[Hx,f1,f2] = freqz2(hx);
figure();
surf(f1,f2,abs(Hx))
axis('equal')
xlabel('x')
ylabel('y')
title('Surface plot of Hx')

[Hy,f1,f2] = freqz2(hy);
figure
surf(f1,f2,abs(Hy))
axis('equal')
xlabel('x')
ylabel('y')
title('Surface plot of Hy')

% b) is in a separate file

% c)
circuit = imread('circuit.tif');
imagedouble = im2double(circuit);
figure();
image(circuit)
title('Original image')
edgeB =imtool(edgeDetect(imagedouble,0.07,2));
set(edgeB,'Name','L2, 0.07');
med = median(imagedouble(:));
edgeM = imtool(edgeDetect(imagedouble,med,2));
set(edgeM,'Name','L2, Median');
%The median does not work, we cannot see anything

edge_1 = imtool(edgeDetect(imagedouble,0.01,2));
edge_2 = imtool(edgeDetect(imagedouble,0.1,2));
set(edge_1,'Name','L2, 0.01');
set(edge_2,'Name','L2, 0.1');
%Observing the behavior

% d)
edgeB2 = imtool(edgeDetect(imagedouble,0.07,1));
set(edgeB2,'Name','L1, 0.07');
edgeM2 = imtool(edgeDetect(imagedouble,med,1));
set(edgeM2,'Name','L1, Median');
edge_3 = imtool(edgeDetect(imagedouble,0.01,1));
set(edge_3,'Name','L1, 0.01');
edge_4 =imtool(edgeDetect(imagedouble,0.1,1));
set(edge_4,'Name','L1, 0.1');