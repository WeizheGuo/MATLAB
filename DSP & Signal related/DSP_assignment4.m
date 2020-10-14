% Weizhe Guo
% ECE 310
% Problem Set 4

%% Problem 2
% a)
fs = 50e3;
N = 256;
bin_spacing = fs/N;

% b)
f = 10e3;
k1 = f/bin_spacing;
k2 = (fs - f)/bin_spacing;

% c)
w_diff = abs(2*pi*(f-bin_spacing*51)/fs);
dir_dc = diric(0, 250);
dir_w = diric(w_diff, 250);
strad_loss = 20*log10(abs(dir_dc/ dir_w));

% d)
ham = hamming(250);
ham_freq = freqz(ham, 1, [0 w_diff]);
ham_strad_loss = 20*log10(abs(ham_freq(1)/abs(ham_freq(2)))); 

% e) 
N_new = 512;
bin_new = fs/N_new;

w_new = abs((2*pi*(f- bin_new*51))/fs);
dir_w_new = diric(w_new, 250);
strad_loss_new = 20*log10(abs(dir_dc/dir_w_new));

ham_freq_new = freqz(ham, 1, [0 w_new]);
ham_strad_loss_new = 20*log10(abs(ham_freq_new(1)/abs(ham_freq_new(2)))); 

%% Problem 6
N = 1000;
n = 1024;
f = 20e6;
fs = 100e6;
r = 30;
t = (0:999)*1/fs;
y = 2*sin(2*pi*f*t);
AWGN = sqrt(0.2)*randn(size(t));
y2 = y + AWGN;
chebw = chebwin(N, r);
win_out = y2 .* transpose(chebw);
f = (-n/2:(n/2-1))*fs/n;
fft_result = fft(win_out, n);
fft_y = 20*log10(abs(fftshift(fft_result)));

figure;
plot(f, fft_y);
xlim([-fs/2 fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Sinewave FFT');

%% Problem 7
% a)
w0 = zeros(1,1000);
for x = 0:249
    w0(4*x+3) = 1;
end
xhat = w0 .* y2;

fft_result2 = fft(xhat, n);
fft_y2 = 20*log10(abs(fftshift(fft_result2)));

figure;
plot(f, fft_y2);
xlim([-fs/2 fs/2]);
xlabel('Hertz (Hz)');
ylabel('Magnitude (dB)');
title('xhat[n] FTT');

% Peaks are at the following freq: -45, -30, -20, -5, 5, 20, 30, 45 (MHz).
% They appear in pairs. 

% b)
k2 = (0:999);
fft_y3 = fft(w0, N);

figure;
stem(k2, abs(fft_y3));
xlabel('Index k');
ylabel('Magnitude');
title('W0[k] Stem');