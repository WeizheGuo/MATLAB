% Weizhe Guo
% ECE 310
% Problem Set 8

%% 4
% a)
rp = 1.5;
rs = 30;
[b_e, a_e] = ellip(4, rp,rs,[0.3 0.6]);
H_z = tf(b_e,a_e,-1);
[h_e,w_e] = freqz(b_e,a_e, 800);

figure
plot(w_e, 20*log10(abs(h_e)));
title('Elliptic filter Magnitude response');
xlabel('Freq');
ylabel('Magnitude (dB)');

% b)
[z_e,p_e,k_e] = tf2zpk(b_e,a_e);
[up_sos, up_g] = zp2sos(z_e,p_e,k_e,'up', 'inf');
[down_sos, down_g] = zp2sos(z_e,p_e,k_e,'down', 'inf');

% c)
[z1,p1,k1] = sos2zp(up_sos(1,:));
[z2,p2,k2] = sos2zp(up_sos(2,:));
[z3,p3,k3] = sos2zp(up_sos(3,:));
[z4,p4,k4] = sos2zp(up_sos(4,:));
up_zeros = [abs(z1) abs(z2) abs(z3) abs(z4)];
up_poles = [abs(p1) abs(p2) abs(p3) abs(p4)];

% the magnitudes of zeros are all 1 and the mags of poles are in increasing
% order and agrees with our expectation

[z5,p5,k5] = sos2zp(down_sos(1,:));
[z6,p6,k6] = sos2zp(down_sos(2,:));
[z7,p7,k7] = sos2zp(down_sos(3,:));
[z8,p8,k8] = sos2zp(down_sos(4,:));
down_zeros = [abs(z5) abs(z6) abs(z7) abs(z8)];
down_poles = [abs(p5) abs(p6) abs(p7) abs(p8)];

% Similarly, here we notice that the mags of zeros are 1 and the mags of
% zeros are in decreasing order

% d)
w=linspace(0,pi,1e4);
h1_up = freqz(up_sos(1,1:3),up_sos(1,4:end),w)*up_g;
h2_up = freqz(up_sos(2,1:3),up_sos(2,4:end),w).*h1_up;
h3_up = freqz(up_sos(3,1:3),up_sos(3,4:end),w).*h2_up;
h4_up = freqz(up_sos(4,1:3),up_sos(4,4:end),w).*h3_up;
figure;
hold on
plot(w/pi,20*log10(abs(h1_up)))
hold on
plot(w/pi,20*log10(abs(h2_up)))
hold on
plot(w/pi,20*log10(abs(h3_up)))
hold on
plot(w/pi,20*log10(abs(h4_up)))
title('Magnitude Response (up)')
legend('H1','H1*H2','H1*H2*H3','H1*H2*H3*H4')
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')

h1_down = freqz(down_sos(1,1:3),down_sos(1,4:end),w)*down_g;
h2_down = freqz(down_sos(2,1:3),down_sos(2,4:end),w).*h1_down;
h3_down = freqz(down_sos(3,1:3),down_sos(3,4:end),w).*h2_down;
h4_down = freqz(down_sos(4,1:3),down_sos(4,4:end),w).*h3_down;
figure;
hold on
grid on
plot(w/pi,20*log10(abs(h1_down)))
hold on
plot(w/pi,20*log10(abs(h2_down)))
hold on
plot(w/pi,20*log10(abs(h3_down)))
hold on
plot(w/pi,20*log10(abs(h4_down)))
title('Magnitude Response (down)')
legend('H1','H1*H2','H1*H2*H3','H1*H2*H3*H4')
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')

% e)
B = [sos_up(1,1:3)./sos_down(4,1:3),sos_up(2,1:3)./sos_down(3,1:3)];
A = [sos_up(1,4:end)./sos_down(4,4:end),sos_up(2,4:end)./sos_down(3,4:end)];

% The denominator equals to each other, the numerator factors in each SOS
% are only off by a scaling factor

%% 5
% a)
b_full = [0.1336, 0.0568, 0.0563, 0.1336];
a_full = [1, -1.5055, 1.2630, -0.3778];

b_1 = [-0.4954, 1];                          
a_1 = [1, -0.4954];

b_2 = [0.7632, -1.0101, 1];
a_2 = [1, -1.0101, 0.7632];

b_fullq = fi(b_full,1,5,3);
b_fullq = b_fullq.data;
a_fullq = fi(a_full,1,5,3);
a_fullq = a_fullq.data;

b_1_q = fi(b_1,1,5,3);
b_1_q = b_1_q.data;
a_1_q = fi(a_1,1,5,3);
a_1_q = a_1_q.data;
b_2_q = fi(b_2,1,5,3);
b_2_q = b_2_q.data;
a_2_q = fi(a_2,1,5,3);
a_2_q = a_2_q.data;

% b)
mask = [1,-1,1,-1];

H_full_dc = sum(b_full)/sum(a_full);
H_full_pi = (sum(b_full.*mask))/(sum(a_full.*mask));

HA_dc = 0.5*(sum(b_1)/sum(a_1) + sum(b_2)/sum(a_2));
HA_pi = 0.5*(sum(b_1.*mask(1:2))/sum(a_1.*mask(1:2)) + sum(b_2.*mask(1:3))/sum(a_2.*mask(1:3)));
H_fullq_dc = sum(b_fullq)/sum(a_fullq);
H_fullq_pi = (sum(b_fullq.*mask))/(sum(a_fullq.*mask));

HA_q_dc = 0.5*(sum(b_1_q)/sum(a_1_q) + sum(b_2_q)/sum(a_2_q));
HA_q_pi = 0.5*(sum(b_1_q.*mask(1:2))/sum(a_1_q.*mask(1:2)) + sum(b_2_q.*mask(1:3))/sum(a_2_q.*mask(1:3)));

H_fullq_dc_error = abs(20*log10(H_fullq_dc/H_full_dc));
HAq_dc_error = abs(20*log10(HA_q_dc/H_full_dc));
H_fullq_pi_error = abs(20*log10(H_fullq_pi/H_full_pi));
HAq_pi_error = abs(20*log10(HA_q_pi/H_full_pi));


% c)
[H_q_a1,w] = freqz(b_1_q,a_1_q,1e4);
H_q_a2 = freqz(b_2_q,a_2_q,w);
H_full = 20*log10(abs(freqz(b_full,a_full,w)));
H_q_full0 = 20*log10(abs(freqz(b_fullq,a_fullq,w)));
H_q_all = 20*log10(abs(0.5*(H_q_a1 + H_q_a2)));

figure
plot(w,H_full);
hold on
plot(w,H_q_full0);
hold on
plot(w,H_q_all);
title('Comparison among three filters')
xlabel('Frequency (rad/sample)')
ylabel('Magnitude (dB)')
legend('Original','Quantized Original','Quantized Allpass')
ylim([-40,0])
xlim([0 pi])


error_q_full = max(abs(H_full-H_q_full0));
error_q_a = max(abs(H_full-H_q_all));

% d)
max_dev_q_0 = max(abs(H_full(1:3001)-H_q_full0(1:3001)));
max_dev_q_a = max(abs(H_full(1:3001)-H_q_all(1:3001)));

ripple_p_0 = max(H_full(1:3001)) - min(H_full(1:3001));
ripple_p_q_0 = max(H_q_full0(1:3001))-min(H_q_full0(1:3001));
ripple_p_q_a = max(H_q_all(1:3001))-min(H_q_all(1:3001));
ripple_s_0 = max(H_full(4070:end))- min(H_full(4103:end));
ripple_s_q_0 = abs(max(H_q_full0(4070:end)))- min(H_q_full0(4070:end));
ripple_s_q_a = abs(max(H_q_all(4070:end)))- min(H_q_all(4070:end));

max_gain_0 = max(H_full(4070:end));
max_gain_q_0 = max(H_q_full0(4070:end));
max_gain_all = max(H_q_all(4070:end));

% q_a has very small deviation from the inf. precision filter in passband
% q_0 has larger deviation 

% The originaal and q_a are close and have a smaller ripple in passband
% The quantized original has a much larger ripple in passband

% For the gain, the original and q_a are close and both satisfy the specs,
% whereas q0 has a larger gain and does not satisfy the spec.

% e)
wp = w(1:3001);
a_qA = conv(a_1_q,a_2_q);
b_qA = conv(b_1_q,a_2_q)+conv(b_2_q,a_1_q);

gpdelay_orig = grpdelay(b_full,a_full,wp);
gpdelay_quant = grpdelay(b_fullq,a_fullq,wp);
gpdelay_all = grpdelay(b_qA,a_qA,wp);

figure
plot(wp,gpdelay_orig);
hold on
plot(wp,gpdelay_quant);
hold on
plot(wp,gpdelay_all);
title('Group Delay of three filters')
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
legend('Original','Quantized Original','Quantized Parallel Config')

% The group delays of allpass filters and quantized full filters are
% identical and smoother than the original filter

% f)
figure()
[z_orig,p_orig,k_orig] = tf2zp(b_full,a_full);
zplane(z_orig, p_orig); 
title('Original Filter');

figure()
[z_quant,p_quant,k_quant] = tf2zp(b_fullq,a_fullq);
zplane(z_quant, p_quant); 
title('Quantized Filter');

figure()
[z_all,p_all,k_all] = tf2zp(b_qA,a_qA);
zplane(z_all, p_all);
title('Quantized Sum of Allpasses')

% The zeros are on the unit circle except for the orginal one.
% No poles are outside of the unit circle.