% Weizhe Guo
% ECE 310
% Problem Set 3

% Analog filter via bilinear transform
wp = [100e3 130e3]*2*pi;
ws = [90e3 140e3]*2*pi;
rp = 2;
rs = 30;
[nb_ana,Wn_ana] = buttord(wp,ws,rp,rs,'s');
[zb_ana,pb_ana,kb_ana] = butter(nb_ana,Wn_ana,'bandpass','s');
[nc1_ana,Wp1_ana] = cheb1ord(wp,ws,rp,rs,'s');
[zc1_ana,pc1_ana,kc1_ana] = cheby1(nc1_ana,rp,Wp1_ana,'bandpass','s');
[nc2_ana,Wp2_ana] = cheb2ord(wp,ws,rp,rs,'s');
[zc2_ana,pc2_ana,kc2_ana] = cheby2(nc2_ana,rs,Wp2_ana,'bandpass','s');
[ne_ana,We_ana] = ellipord(wp,ws,rp,rs,'s');
[ze_ana,pe_ana,ke_ana] = ellip(ne_ana,rp,rs,We_ana,'bandpass','s');

% Degital via bilinear transform
fs = 400e3;
wp_dig = [100e3 130e3]*2/fs; % normalization
ws_dig = [90e3 140e3]*2/fs;
[nb_dig,Wn_dig] = buttord(wp_dig,ws_dig,rp,rs);
[zb_dig,pb_dig,kb_dig] = butter(nb_dig,Wn_dig,'bandpass');
[nc1_dig,Wp1_dig] = cheb1ord(wp_dig,ws_dig,rp,rs);
[zc1_dig,pc1_dig,kc1_dig] = cheby1(nc1_dig,rp,Wp1_dig,'bandpass');
[nc2_dig,Wp2_dig] = cheb2ord(wp_dig,ws_dig,rp,rs);
[zc2_dig,pc2_dig,kc2_dig] = cheby2(nc2_dig,rs,Wp2_dig,'bandpass');
[ne_dig,We_dig] = ellipord(wp_dig,ws_dig,rp,rs);
[ze_dig,pe_dig,ke_dig] = ellip(ne_dig,rp,rs,We_dig,'bandpass');

% Digital via impulse invariance
[bb_ana, ab_ana] = butter(nb_ana,Wn_ana,'s');
[bc1_ana, ac1_ana] = cheby1(nc1_ana,rp,Wp1_ana,'s');
[bc2_ana, ac2_ana] = cheby2(nc2_ana,rs,Wp2_ana,'s');
[be_ana,ae_ana] = ellip(ne_ana,rp,rs,We_ana,'s');

[bb_dig,ab_dig] = impinvar(bb_ana,ab_ana,fs);
[bc1_dig,ac1_dig] = impinvar(bc1_ana,ac1_ana,fs);
[bc2_dig,ac2_dig] = impinvar(bc2_ana,ac2_ana,fs);
[be_dig,ae_dig] = impinvar(be_ana,ae_ana,fs);

% Zplane plot
figure
subplot(2,2,1);
zplane(zb_ana,pb_ana)
grid
title('Butterworth Filter (Analog)')
subplot(2,2,2);
zplane(zc1_ana,pc1_ana)
grid
title('Cheb1 Filter (Analog)')
subplot(2,2,3);
zplane(zc2_ana,pc2_ana)
grid
title('Cheb2 Filter (Analog)')
subplot(2,2,4);
zplane(ze_ana,pe_ana)
grid
title('Elliptic Filter (Analog)')
hold on;

figure
subplot(2,2,1);
zplane(zb_dig,pb_dig)
grid
title('Butterworth Filter (Digital)')
subplot(2,2,2);
zplane(zc1_dig,pc1_dig)
grid
title('Cheb1 Filter (Digital)')
subplot(2,2,3);
zplane(zc2_dig,pc2_dig)
grid
title('Cheb2 Filter (Digital)')
subplot(2,2,4);
zplane(ze_dig,pe_dig)
grid
title('Elliptic Filter (Digital)')
hold on;

figure
subplot(2,2,1);
zplane(bb_dig,ab_dig)
grid
title('Butterworth Filter (Digital via impulse inv)')
subplot(2,2,2);
zplane(bc1_dig,ac1_dig)
grid
title('Cheb1 Filter (Digital via impulse inv)')
subplot(2,2,3);
zplane(bc2_dig,ac2_dig)
grid
title('Cheb2 Filter (Digital via impulse inv)')
subplot(2,2,4);
zplane(be_dig,ae_dig)
grid
title('Elliptic Filter (Digital via impulse inv)')
hold on;

% Magnitude response plot
% Analog
w = linspace(0,200e3*2*pi,2000);
hb_ana = freqs(bb_ana,ab_ana,w);
figure;
subplot(2,2,1);
plot(w,20*log10(abs(hb_ana)));
axis([0 200e3*2*pi -50 1]);
yticks(-50:5:1);
title('Butterworth Magnitude (Analog)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

hc1_ana = freqs(bc1_ana,ac1_ana,w);
subplot(2,2,2);
plot(w,20*log10(abs(hc1_ana)));
axis([0 200e3*2*pi -50 1]);
yticks(-50:5:1);
title('Chebyshev I Magnitude (Analog)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

hc2_ana = freqs(bc2_ana,ac2_ana,w);
subplot(2,2,3);
plot(w,20*log10(abs(hc2_ana)));
axis([0 200e3*2*pi -50 1]);
yticks(-50:5:1);
title('Chebyshev II Magnitude (Analog)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

he_ana = freqs(be_ana,ae_ana,w);
subplot(2,2,4);
plot(w,20*log10(abs(he_ana)));
axis([0 200e3*2*pi -50 1]);
yticks(-50:5:1);
title('Elliptic Magnitude (Analog)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

% Digital via impulse invariance
f = linspace(0, fs/2, 2000);
hb_dig = freqz(bb_dig,ab_dig,f,fs);
figure;
subplot(2,2,1);
plot(f,20*log10(abs(hb_dig)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Butterworth Magnitude (Digital via impulse invariance)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

hc1_dig = freqz(bc1_dig,ac1_dig,f,fs);
subplot(2,2,2);
plot(f,20*log10(abs(hc1_dig)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Chebyshev I Magnitude (Digital via impulse invariance)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

hc2_dig = freqz(bc2_dig,ac2_dig,f,fs);
subplot(2,2,3);
plot(f,20*log10(abs(hc2_dig)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Chebyshev II Magnitude (Digital via impulse invariance)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

he_dig = freqz(be_dig,ae_dig,f,fs);
subplot(2,2,4);
plot(f,20*log10(abs(he_dig)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Elliptic Magnitude (Digital via impulse invariance)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

% Digital
[bb_dig2, ab_dig2] = butter(nb_dig,Wn_dig);
[bc1_dig2, ac1_dig2] = cheby1(nc1_dig,rp,Wp1_dig);
[bc2_dig2, ac2_dig2] = cheby2(nc2_dig,rs,Wp2_dig);
[be_dig2,ae_dig2] = ellip(ne_dig,rp,rs,We_dig);

hb_dig2 = freqz(bb_dig2,ab_dig2,f,fs);
figure;
subplot(2,2,1);
plot(f,20*log10(abs(hb_dig2)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Butterworth Magnitude (Digital)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

hc1_dig2 = freqz(bc1_dig2,ac1_dig2,f,fs);
subplot(2,2,2);
plot(f,20*log10(abs(hc1_dig2)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Chebyshev I Magnitude (Digital)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

hc2_dig2 = freqz(bc2_dig2,ac2_dig2,f,fs);
subplot(2,2,3);
plot(f,20*log10(abs(hc2_dig2)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Chebyshev II Magnitude (Digital)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

he_dig2 = freqz(be_dig2,ae_dig2,f,fs);
subplot(2,2,4);
plot(f,20*log10(abs(he_dig2)));
axis([0 fs/2 -50 1]);
yticks(-50:5:1);
title('Elliptic Magnitude (Digital)');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

% Report actual order
% actual order = 2 * returned order from matlab function
fprintf('Butterworth filter have different orders in two cases: Digital Butterworth has order 16, whereas analog case has order 18\n')
fprintf('Both Cheb1 filter and Cheb2 filter have the same order -- 10 for two cases\n')
fprintf('Ellip filter has the same order -- 6 for two cases\n')

% Calculate time constant for each system
% Find dominant pole and than find the time const of the pole
% Analog
tb_ana = 1/abs(max(real(pb_ana)));
tc1_ana = 1/abs(max(real(pc1_ana)));
tc2_ana = 1/abs(max(real(pc2_ana)));
te_ana = 1/abs(max(real(pe_ana)));
% Digital
T = 1/fs;
tb_dig = T/abs(log(max(abs(pb_dig))));
tc1_dig = T/abs(log(max(abs(pc1_dig))));
tc2_dig = T/abs(log(max(abs(pc2_dig))));
te_dig = T/abs(log(max(abs(pe_dig))));
% Digital via impulse
[zb2_dig,pb2_dig,kb2_dig] = tf2zp(bb_dig,ab_dig);
[zc12_dig,pc12_dig,kc12_dig] = tf2zp(bc1_dig,ac1_dig);
[zc22_dig,pc22_dig,kc22_dig] = tf2zp(bc2_dig,ac2_dig);
[ze2_dig,pe2_dig,ke2_dig] = tf2zp(be_dig,ae_dig);

tb2_dig = T/abs(log(max(abs(pb2_dig))));
tc12_dig = T/abs(log(max(abs(pc12_dig))));
tc22_dig = T/abs(log(max(abs(pc22_dig))));
te2_dig = T/abs(log(max(abs(pe2_dig))));
% As we see above, the impulse invariance perserves time constant better
% than bilinear transform

% Do bilinear by hand
wp_dig2 = [100e3 130e3]*2*pi/fs; 
ws_dig2 = [90e3 140e3]*2*pi/fs;
omega_p = tan(wp_dig2/2);
omega_s = tan(ws_dig2/2);
omega0 = sqrt(omega_p(1)*omega_p(2));
B = omega_p(1) - omega_p(2);
omega_stop = min([(abs((omega_s(1)^2-omega0^2)/(B*omega0))) (abs((omega_s(2)^2-omega0^2)/(B*omega0)))]);
