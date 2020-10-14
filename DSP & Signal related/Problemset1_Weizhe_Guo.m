% Weizhe Guo
% ECE-211
% 02/07/2018
% Problem Set 1

%% Problem 5
% d)
t = linspace(0,1,1001);
approxS = 0;
for m = -5:5
    approxS = approxS+1/(2*pi*(1+1j*m))*(1-exp(-2*pi*(1+1j*m)))*exp(1j*m*2*pi*t);
end
maximagS = 1j*max(abs(imag(approxS)))
realS = real (approxS);
S = exp(-2*pi*t);

figure
plot(t,realS);
xlabel('time');
ylabel('signal');
title('comparison between original and approximate signal');
grid on;
hold on;

plot(t,S);
legend('approximate signal','original signal');

% e)
n = linspace(-5,5,11);
Cn = (1-exp(-2*pi))*(1-1j*n)./(2*pi*(n.^2+1));
figure
stem(n,real(Cn))
xlabel('n(discrete)');
ylabel('fourier coefficient');
title('fourier coefficient for different n');
hold on;