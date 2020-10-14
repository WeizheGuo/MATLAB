% Weizhe Guo
% ECE 310
% Problem Set 2

% a) Calculate w0 and B
wphi = 130000*2*pi;
wplo = 100000*2*pi;
wslo = 90000*2*pi;
wshi = 140000*2*pi;
B = wphi - wplo;
w0 = sqrt(wphi*wplo);

%select the lower one after mapping
wproto_lo = abs((wslo^2-w0^2)/(B*wslo));
wproto_hi = abs((wshi^2-w0^2)/(B*wshi));
%select wproto_hi as it is smaller

% b) calculate n
n_butter = 0.5*log10((10^(30/10)-1)/(10^(2/10)-1))/log10(wproto_hi/1);
n_cheby = acosh(sqrt((10^(30/10)-1)/(10^(2/10)-1)))/acosh(wproto_hi/1);
% The order should be the integer greater than the above calculated values
% n_butter = 9 and n_cheby = 5

% c) use matlab built-in function
wp = [100e3 130e3];
ws = [90e3 140e3];
rp = 2;
rs = 30;
[nb,Wn] = buttord(wp,ws,rp,rs,'s');
[zb,pb,kb] = butter(nb,Wn,'s');
[nc1,Wp1] = cheb1ord(wp,ws,rp,rs,'s');
[zc1,pc1,kc1] = cheby1(nc1,rp,Wp1,'s');
[nc2,Wp2] = cheb2ord(wp,ws,rp,rs,'s');
[zc2,pc2,kc2] = cheby2(nc2,rs,Wp2,'s');
[ne,We] = ellipord(wp,ws,rp,rs,'s');
[ze,pe,ke] = ellip(ne,rp,rs,We,'s');

% d) zplane plot
figure
subplot(2,2,1);
zplane(zb,pb)
grid
title('Butterworth Filter')
subplot(2,2,2);
zplane(zc1,pc1)
grid
title('Cheb1 Filter')
subplot(2,2,3);
zplane(zc2,pc2)
grid
title('Cheb2 Filter')
subplot(2,2,4);
zplane(ze,pe)
grid
title('Elliptic Filter')
hold on;

% e) prototype lowpass filter
figure;
[ne1,We1] = ellipord(1,wproto_hi,rp,rs,'s');
[ze1,pe1,ke1] = ellip(ne1,rp,rs,wp,'s');
zplane(ze1,pe1);
grid
title('Elliptic Lowpass Filter')
hold on;

% f)
w = linspace(0,200e3*2*pi,2000);
[bb, ab] = butter(nb,Wn,'s');
hb = freqs(bb,ab,w);
figure;
subplot(2,1,1);
plot(w,db(hb));
axis([0 200e3 -100 20]);
yticks(-100:10:20);
title('Butterworth Magnitude');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
hline = refline([0 -2]);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline([0 -30]);
hline.Color = 'b';
hline.LineStyle = '--';

subplot(2,1,2);
plot(w,unwrap(angle(hb))*180/pi);
title('Butterworth Phase');
ylabel('Phase (degree)');
xlabel('Frequency (Hz)');

[bc1, ac1] = cheby1(nc1,rp,Wp1,'s');
hc1 = freqs(bc1,ac1,w);
figure;
subplot(2,1,1);
plot(w,db(hc1));
axis([0 200e3 -100 20]);
yticks(-100:10:20);
title('Chebyshev I Magnitude');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
hline = refline([0 -2]);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline([0 -30]);
hline.Color = 'b';
hline.LineStyle = '--';

subplot(2,1,2);
plot(w,unwrap(angle(hc1))*180/pi);
title('Chebyshev I Phase');
ylabel('Phase (degree)')
xlabel('Frequency (Hz)');

[bc2, ac2] = cheby2(nc2,rs,Wp2,'s');
hc2 = freqs(bc2,ac2,w);
figure;
subplot(2,1,1);
plot(w,db(hc2));
axis([0 200e3 -100 20]);
yticks(-100:10:20);
title('Chebyshev II Magnitude');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
hline = refline([0 -2]);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline([0 -30]);
hline.Color = 'b';
hline.LineStyle = '--';

subplot(2,1,2);
plot(w,unwrap(angle(hc2))*180/pi);
title('Chebyshev II Phase');
ylabel('Phase (degree)');
xlabel('Frequency (Hz)');

[be,ae] = ellip(ne,rp,rs,We,'s');
he = freqs(be,ae,w);
figure;
subplot(2,1,1);
plot(w,db(he));
axis([0 200e3 -100 20]);
yticks(-100:10:20);
title('Elliptic Magnitude');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
hline = refline([0 -2]);
hline.Color = 'b';
hline.LineStyle = '--';
hline = refline([0 -30]);
hline.Color = 'b';
hline.LineStyle = '--';

subplot(2,1,2);
plot(w,unwrap(angle(he))*180/pi);
title('Elliptic Phase');
ylabel('Phase (degree)');
xlabel('Frequency (Hz)');

% g) Superimpose and zoom in
figure();
plot(w,db(hb));
hold on;
plot(w,db(hc1));
hold on;
plot(w,db(hc2));
hold on;
plot(w, db(he));
hold on;
xlim([100000 130000]);
ylim([-2 0]); %zoom into passband
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');

% h) Attenuation at stopband edge
ws_rad = [90e3 140e3];
hb_stop = db(freqs(bb,ab,ws_rad));
hc1_stop = db(freqs(bc1,ac1,ws_rad));
hc2_stop = db(freqs(bc2,ac2,ws_rad));
he_stop = db(freqs(be,ae,ws_rad));
Filter = {'ButterWorth';'Chebyshev I';'Chebyshev II';'Elliptic'};
Attenuation = [hb_stop; hc1_stop; hc2_stop; he_stop];
Table = table(Filter,Attenuation);