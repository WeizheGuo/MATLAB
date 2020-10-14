% Weizhe Guo & Zheng Liu
% ECE302-Bonus

clear all;
clc;
load pj2data.mat

% I did several questions until I realize I dont have to do them

%% Part A: Autocorrelation and MATLAB
y1 = y(1:32);
n = linspace(1,63,63);
cyy = xcorr(y1,y1,'biased');

figure;
plot(n, cyy);
xlabel('n')
ylabel('Autocorrelation of y1[n]')
title('Autocorrelation of y1[n] Using xcorr (biased)')

cyy_conv = conv(y1,fliplr(y1));
figure;
plot(n,cyy_conv);
xlabel('n')
ylabel('Autocorrelation of y1[n]')
title('Autocorrelations of y1[n] Using conv')

% Part A.1)
cyy_xcorr_unbiased = xcorr(y1,y1,'unbiased');
figure;
plot(n, cyy_xcorr_unbiased); % plot the unbiased result
xlabel('n')
ylabel('Autocorrelation of y1[n] (unbiased)')
title('Autocorrelation of y1[n] Using xcorr (unbiased)')
% The unbiased result will reduce the denominator

% Part A.2.a)
% The Fourier Transform of the Autocorrelation Function is the Power
% Spectrum, which is non-negative and real. 

%% C
% C.1
cyy_corr = xcorr(y,y,'biased');
[coef_a_2, est_err_2, PDS_est_2] = Calculate_PDS_est(cyy_corr, 2);
[coef_a_3, est_err_3, PDS_est_3] = Calculate_PDS_est(cyy_corr, 3);
[coef_a_4, est_err_4, PDS_est_4] = Calculate_PDS_est(cyy_corr, 4);
[coef_a_5, est_err_5, PDS_est_5] = Calculate_PDS_est(cyy_corr, 5);
[coef_a_6, est_err_6, PDS_est_6] = Calculate_PDS_est(cyy_corr, 6);
[coef_a_7, est_err_7, PDS_est_7] = Calculate_PDS_est(cyy_corr, 7);

%[freq_resp, w] = freqz(ones(2,1), coef_a_2, 512);
%PDS_est = abs(freq_resp).^2;

figure;
n = linspace(1,512,512);
plot(n,PDS_est_2);
hold on
plot(n,PDS_est_3);
plot(n,PDS_est_4);
plot(n,PDS_est_5);
plot(n,PDS_est_6);
plot(n,PDS_est_7);
xlabel('n')
ylabel('PDS estimation')
title('PDS estimation made for different orders')
legend('2nd order','3rd order','4th order', '5th order', '6th order', '7th order') 


function [coef_a, est_err, PDS_est] = Calculate_PDS_est(cyy_corr, ord_num) 
    [coef_a, est_err] = levinson(cyy_corr, ord_num);
    % calculate freq response
    [freq_resp, w] = freqz(ones(ord_num,1), coef_a, 512);
    % square of magnitude of freq response is the estimation
    PDS_est = abs(freq_resp).^2;
end