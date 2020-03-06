clc
clear all
% close all

% Htop = 3*Ht
Data = ;
d = Data(:,1);neff = (Data(:,2));perm = neff.^2;

figure;
yyaxis right;
plot(d*1e9,abs(real(perm)));
set(gcf, 'Position', [00, 00, 400, 300]);
ylabel('Real permittivity');
xlabel('Waveguide width');
yyaxis left;
hold on;
plot(d*1e9,imag(perm));
ylabel('Imaginary permittivity');
xlabel('Waveguide width');
legend(["Imaginary permittivity","Real permittivity"]);
title('PECVD 50nm SiN');
xlim([150 500])

k0 = 2*pi/1550e-5;
L = 20e-6;
M = 300;
ne = 1.8;
z = -L/2:(L/(M - 1)):L/2;
wz = 2*ne*z/L;
permz = wz.^2;
figure;
plot(z,permz);

for ind = 1:1:M
    [~,inx] = min(abs((real(permz(ind) - perm)))); 
    d_Al(ind) = d(inx); 
    eps(ind) = perm(inx); 
end

figure;
yyaxis left;
plot(z*1e6,real(eps));
set(gcf, 'Position', [00, 00, 400, 300]);
xlabel('x(\mum)');
ylabel('Re[\epsilon]');
yyaxis right
plot(z*1e6,imag(eps));
set(gcf, 'Position', [00, 00, 400, 300]);
xlabel('x(\mum)');
ylabel('Im[\epsilon]');

figure;
plot(z*1e6,d_Al*1e9);
set(gcf, 'Position', [00, 00, 400, 300]);
xlabel('x(\mum)');
ylabel('Waveguide width(nm)')

figure;
plot(z,imag(eps))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Comsol Geometry
% a = 0.181/0.161 - 1;
a = 0.161/0.161 - 1;
gc = -a*(abs(2*z/L)).^2 + a + 1;
Al_Y = [d_Al.*gc*1e6 -d_Al.*gc*1e6]/2;
Al_Y_h = Al_Y';
Al_Z = [z -z]*1e6;
Al_Z_h = Al_Z';
% Al_Y = [d_Al-0.0e-6 -d_Al+0.0e-6]/2*1e6;Al_Z = [z -z]*1e6;
% PumpSweep = unique(Lx_Al_Com)*1e9;
SiO_Y = [d_Al + 200e-9 -d_Al - 200e-9]/2*1e6;
SiO_Z = [z -z];


% % %%%%3.3^2 +(0.03*exp(-((z-5.168[um])/2[um])^2) + 0.03*exp(-((z+5.168[um])/2[um])^2) + 0.017*exp(-(z/10[um])^2))*i
% % tm = 0.03*exp(-((z-5.168e-6)/2e-6).^2) + 0.03*exp(-((z+5.168e-6)/2e-6).^2) + 0.017*exp(-(z/10e-6).^2);
% % figure
% % plot(z,tm)