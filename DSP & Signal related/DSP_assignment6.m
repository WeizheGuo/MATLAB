% Weizhe Guo
% ECE 310
% Problem Set 5

%% Problem 1
% a) 
L = 31;
r = 30;
cheb_w = 1/L * chebwin(L,r);
norm_chebw = cheb_w/ sum(cheb_w); % normalization: DC = 1; all coeff sum up to be 1

% b)
figure();
zplane(norm_chebw.');
z0 = roots(norm_chebw.');
w = log(z0)/(1*j);
w0 = min(w(w > 0));
main_width = 2*w0;
rel_width = main_width/(4 * pi / L);

% c)
[cheb_W,f] = freqz(norm_chebw,1,1000,2); 
figure();
plot(f,20*log10(abs(cheb_W)));
title('Magnitude Response');
xlabel('Normalized Frequency');
ylabel('Normalized Magnitude (dB)');

% d)
beta = 3.13;
% kaiser_w = 1/L * kaiser(L, beta);
% norm_kaiserw = kaiser_w/sum(kaiser_w);
% [kaiser_W,f] = freqz(norm_kaiserw,1,1000,2);
% figure();
% plot(f,abs(cheby_W));
% hold on
% plot(f,abs(kaiser_W));
% legend('Cheby Window', 'Kaiser Window');

% e)
figure();
plot(f,20*log10(abs(cheby_W)))
hold on
plot(f,20*log10(abs(kaiser_W)))
title('Magnitude Responses')
legend('Cheb window','Kaiser window')
xlabel('Normalized Frequency')
ylabel('Normalized Magnitude (dB)')

figure();
zplane(norm_chebw.',norm_kaiserw.');
title('Plot of Zeros')
legend('Cheby Zeros','Kaiser Zeros');

% f)
% beta is already reported in d)

% g)
z0_2 = roots(norm_kaiserw.');
w_2 = log(z0_2)/(1*j);
w0_2 = min(w_2(w_2 > 0));
main_width_2 = 2*w0_2;
rel_width_2 = main_width_2/(4 * pi / L);

omega_2 = linspace(w0_2, pi, 1000);
W_max = max(abs(polyval(norm_kaiserw.', exp(1*j*omega_2)))); % find max value of magnitude of W
side_level_2 = min(20*log10(1/W_max));

% h)
energy_total = sum((abs(norm_chebw)).^2);
omega = linspace(real(w0), pi, 1000);
neg_omega = linspace(-real(w0),-pi,1000); 
% w0 is real+ 0.0000i, but matlab thinks it is a complex matrix
% In order to plug in freqz, I make it to be real by doing real(w0)... But actually real(w0) = w0
W_pos = freqz(norm_chebw,1,omega);
W_neg = freqz(norm_chebw,1,neg_omega);
dw = (pi - real(w0))/1000;
energy_side = sum(abs(W_pos).^2*dw+abs(W_neg).^2*dw);
energy_frac = energy_side / energy_total;

energy_total2 = sum((abs(norm_kaiserw)).^2);
neg_omega_2 = linspace(-real(w0_2),-pi,1000);
omega_2 = real(omega_2);
% similar reason to do real()...
W_pos2 = freqz(norm_kaiserw,1,omega_2);
W_neg2 = freqz(norm_kaiserw,1,neg_omega_2);
dw2 = (pi - real(w0_2))/1000;
energy_side2 = sum(abs(W_pos2).^2*dw2+abs(W_neg2).^2*dw2);
energy_frac2 = energy_side2 / energy_total2;

%% Problem 2
% b)
Ap = 2;
As = 30;
delta_p = (10^(Ap/20)-1)/(10^(Ap/20)+1);
delta_s = 10^(-As/20);
dev = [delta_s, delta_p, delta_s];

% c)
f_bound = 1e6*[1.5, 2, 3, 3.5];
f_sample = 10e6;
mag_desired = [0,1,0];

[n,fo,ao,w] = firpmord(f_bound,mag_desired,dev,f_sample);
n = n+4;
pm_filt = firpm(n,fo,ao,w);
% Thanks Alex Zheng Liu for helping me figure out the confusion here. The
% initial design without adding order eventually ends up with a filter not
% satifying the spec. The passband attenuation doesnot satisfy the requirement. 
% Adding 4 to the order was foound to be the optimal design satisfying the 
% spec.
[H_pm,f_pm] = freqz(pm_filt,1,1000, f_sample);

[n2, Wn2, beta2, ftype2] = kaiserord(f_bound, mag_desired, dev, f_sample);
n2 = n2 + rem(n2,2);
k_filt = fir1(n2, Wn2, ftype2, kaiser(n2 + 1, beta2),'noscale');
[H_k,f_k] = freqz(k_filt, 1, 1000, f_sample);
% After I change the order, there is no over-design. PM filter's passband
% is a bit wider

% d)
pm_filt_add0 = [zeros(1,4) pm_filt zeros(1,5)]; 
% the difference between two orders is 9, so adding 4 and 5 columns to the
% beginning and the end
figure();
subplot(2,1,1);
stem(0:max([n n2]),pm_filt_add0);
subplot(2,1,2);
stem(0:max([n n2]),k_filt);


figure();
subplot(2,1,1);
plot(f_pm, 20*log10(abs(H_pm)));
title('PM Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
hold on
ylim([-50,2]);
hline = refline(0,1);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline(0,-1);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline(0, -30);
hline.Color = 'b';
hline.LineStyle = '--';

subplot(2,1,2);
plot(f_k, 20*log10(abs(H_k)));
title('Kaiser Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
hold on
ylim([-50,2]);
hline = refline(0,0);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline(0,-2);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline(0, -30);
hline.Color = 'b';
hline.LineStyle = '--';

%e
result = w.'.*dev;
% Times w and dev together and founds that the product is the same