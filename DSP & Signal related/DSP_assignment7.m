% Weizhe Guo
% ECE 310
% Problem Set 7

%% 4
% The function creates a FIR filter with coeff vector h
syms z
h2FIR = @(h,z) expand(poly2sym(h,z)/z^(length(h)-1));
% The function for finding the paraconjugate is at the bottom

%% 5
% a)
coef0 = [0.15774243 0.69950381 1.06226376 0.44583132 -0.31998660 -0.18351806 0.13788809 0.03892321 -0.04466375 -7.83251152e-4 6.75606236e-3 -1.52353381e-3];
coef1 = fliplr(coef0);
coef3 = coef1;
coef4= coef0;
coef1(2:2:end) = -coef1(2:2:end);
coef4(2:2:end) = -coef4(2:2:end);

w = linspace(0,pi,1000);
H0 = freqz(coef0, 1, w);
H1 = freqz(coef1, 1, w);
figure()
plot(w, abs(H0));
hold on
plot(w, abs(H1)); 
legend("mag of H0(w)", "mag of H1(w)");
title('Magnitude Response of H0 and H1');
xlabel('Frequency(rad)');
ylabel('Magnitude');

% b)
F0 = freqz(coef3, 1, w);
F1 = freqz(coef4, 1, w);
figure()
plot(w, abs(F0));
hold on
plot(w, abs(F1)); 
legend("mag of F0(w)", "mag of F1(w)");
title('Magnitude Response of F0 and F1')
xlabel('Frequency(rad)')
ylabel('Magnitude')
% I think there is a typo:) You want us to compare F0 to H0, F1 to H1
% The two functions have the same shape, and for some selected values, the
% result is almost the same with some small erros

% c)
sum = (abs(H1)).^2 + (abs(H0)).^2;
% The elements of the matrix sum are found to be a constant. The constant
% is 3.994

% d)
H_flip = flip(H0);
figure()
plot(w,abs(H1));
hold on
plot(w,abs(H_flip)); 
legend("H1", "Flipped H0");
% From the graph, we observe that the two curves coincide, and thus proves
% that H0 and H1 are mirror filters

% e)
syms z
E_00 = h2FIR(coef0(1:2:end),z);
E_01 = h2FIR(coef0(2:2:end),z);
E_10 = h2FIR(coef1(1:2:end),z);
E_11 = h2FIR(coef1(2:2:end),z);

E = [E_00,E_01;E_10,E_11];
E_tilde = para(E);
E_prod = simplify(E_tilde*E);
% The result is a diaganol matrix and thus confirms our hypothesis

% f)
R_00 = h2FIR(coef3(2:2:end),z);
R_01 = h2FIR(coef4(2:2:end),z);
R_10 = h2FIR(coef3(1:2:end),z);
R_11 = h2FIR(coef4(1:2:end),z);

R = [R_00,R_01;R_10,R_11];
RE_prod = simplify(R*E);

% function for question 4
function paraconj = para(symFIR)
    syms z
    conju = symFIR';
    paraconj = subs(conju, conj(z) , z^-1);
end