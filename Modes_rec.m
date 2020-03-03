clc
clear all
close all

global L ne
N = 10000;
ne = 3.3; L = 1; k = 10 * pi / L; h = L / N; x = -0.5 * L:h:0.5 * L;
% k = 2*pi/L;
Ksi = diag((2 * ne * x / L).^2 + 2i * 1 * ne / k / L);
ks1 = diag(Ksi);
H0 = diag(2 * ones(N + 1, 1)) - diag(ones(N, 1), -1) - diag(ones(N, 1), 1);
% H0(1,1) = (2 - 0.0i);H0(N+1,N+1) = (2 - 0.0i);
H0(1, 1) = (2 - exp(-1i * ne * k * h));
H0(N + 1, N + 1) = (2 - exp(-1i * ne * k * h));
Hmm = H0 - h * h * k * k * Ksi .* eye(N + 1);
% Hmm(1,1) = real(Hmm(1,1)) +0.01i;Hmm(N+1,N+1) = real(Hmm(N+1,N+1))+0.01i;

[a, b] = eig(Hmm);
figure;
subplot(2, 1, 1);
plot(real(diag(b)), imag(diag(b)), 'b*');
subplot(2, 1, 2);

for kk = 1:1:N + 1
    plot(abs(a(:, kk)));
    hold on;
end

[c, d] = min(std(abs(a)));
figure
subplot(2, 1, 1);
plot(abs(a(:, d)));
hold on;
subplot(2, 1, 2);
plot(angle(a(:, d)));
hold on;
