% Weizhe Guo
% ECE 310
% Problem Set 1

%% Problem 2

% calculate coefficient of H, A, and Hmin
a = [2, 11, 12, -9];
b = [3, 19, 20];
a_A = [27, 189, 297, 139, 20];
b_A = [20, 139, 297, 189, 27, 0];
a_min = [9, 3/2, -2, -1/2];
b_min = [10, 19/2, 3/2];

[phi, w] = phasez(b, a);
[phi_A, w_A] = phasez(b_A, a_A);
[phi_min, w_min] = phasez(b_min, a_min);

figure();
plot (w, phi*pi/180);
hold on;
plot (w_A, phi_A*pi/180);
hold on;
plot (w_min, phi_min*pi/180);
hold on;
legend('H', 'All pass', 'min phase');
xlabel('freq');
ylabel('phase');