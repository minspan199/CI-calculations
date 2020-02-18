clc;
clear all;
close all

N = 5000;
n0 = 3.3;
L = 1.000;
k = 4*pi/L;
lam0 = 1.550;
P = 1;
S0 = 1.000;
Space = P*L/N;
x = -P*L/2:Space:P*L/2;
dd = 50*(rand(1,N+1)-0.5);
% dd = 0
Sx = lam0/2*(n0*n0 - (n0*n0 - lam0*lam0/S0/S0/4)*4*x.*x/P/P/L/L).^(-1/2) + 0.0*(rand(1,N+1) - 0.5);
S0 = Sx(1);
Ky = pi./Sx;
% n1 = 3.287;
% n1 = 3.165;
n1 = n0;
Kis_eff = n1*n1 - lam0*lam0*Ky.*Ky/4/pi/pi;
% Kis_eff = n1*n1*ones(1,N + 1);
ne = sqrt(n1*n1 - lam0*lam0/S0/S0/4);
Xlimit = 190; Ylimt = 700;

% wx = ne*tanh(x*k/2);
% Ksi = (wx).^2 - 1i*ne/2*(sech(x*k/2)).^2;

figure
% plot(x,Kis_eff,'r*')
% hold on
plot(x,Kis_eff,'b--', 'LineWidth',3)
set(gca,'FontSize', 14)
set(gcf, 'Position', [00, 00, 400, 300])

Dth = 2*ne/k/P/L;

parfor X = 1:1:Xlimit
    v22 = zeros(1,Ylimt);
    v12 = zeros(1,Ylimt);
    v21 = zeros(1,Ylimt);
    rPer = 0.2/L;
    for Y = 1:1:Ylimt
        
        Ms = eye(2);
        Km = (0.01*X + 0.001i*(Y-Ylimt/2))*k;
        Kis_imag = (rPer*rPer*Dth/((k - (Km))^2 + rPer*rPer));
%         Kis_imag = Dth;
        Ksi = Kis_eff - 1i*Kis_imag;
        n = sqrt(Ksi);
        rPer = 0.2/L;
        D0 = [1 1;ne*Km -ne*Km];
        Dn = [1 1;ne*Km -ne*Km];
        for kt = 1:1:N+1
            Dv = [1 1;n(kt)*Km -n(kt)*Km];
            Pv = [exp(1i*n(kt)*Km*Space) 0;0 exp(-1i*n(kt)*Km*Space)];
            Mv = Dv*Pv*Dv^(-1);
            Ms = Ms*Mv;
        end
        
        Ms = D0^(-1)*Ms*Dn;
        v22(1,Y) = abs(Ms(2,2));
        v12(1,Y) = abs(Ms(1,2));
        v21(1,Y) = abs(Ms(2,1));
    end
    Z_22(X,:) = v22;
    Z_12(X,:) = v12;
    Z_21(X,:) = v21;
    X
end

figure
A = flip(log(abs(1./Z_22')));
imagesc(0.01*(1:1:Xlimit),0.001*((Ylimt:-1:1) - Ylimt/2), A, [-4 7])
set(gca,'YDir','normal')
set(gcf, 'Position', [00, 00, 400, 300])
set(gca,'FontSize', 14)
colorbar
%% Data Processing


Dth = 2*ne/k/P/L;
Km = (1 - 1i*0.0)*k;
Ksi = Kis_eff - 1i*Dth;%%*(sech(x)-0.5)
n = sqrt((Ksi));
% Km = (0.57 - 1i*0.013)*k;
% Km = (1.38 + 1i*0.12)*k;
Ms = eye(2);
D0 = [1 1;ne*Km -ne*Km];
Dn = [1 1;ne*Km -ne*Km];
Dv = [1 1;n(N+1)*Km -n(N+1)*Km];
out(:,1) = [0;1];
out(:,2) = Dv^(-1)*Dn*out(:,1);
for kt = N-1:-1:2
    Dv = [1 1;n(kt) -n(kt)];
    Dv_ = [1 1;n(kt-1) -n(kt-1)];
    Pv = [exp(1i*n(kt)*Space*Km) 0;0 exp(-1i*n(kt)*Space*Km)];
    Mv_ = Dv_^(-1)*Dv*Pv;
    out(:,N - kt + 2) = Mv_*out(:,N - kt + 1);
end


figure
in = 1;
for kt = 1:in:N
    Psi((kt + in - 1)/in) = out(1,kt)*exp(-1i*n(kt)*Km*Space/2) + ...
        out(2,kt)*exp(1i*n(kt)*Km*x(kt)*Space/2);
end
plot(abs(Psi),'b*')
set(gca,'FontSize', 14)
set(gcf, 'Position', [00, 00, 400, 300])
set(gca,'xtick',[0,N/4,N/2,N/4*3,N]);
set(gca,'xticklabel',{'-0.5','-0.25','0','0.25','0.5'});
axis([-0 N 0 1.2])

figure
plot(abs(out(1,:)),'b*')
hold on
plot(abs(out(2,:)),'r*')
set(gca,'FontSize', 14)
set(gcf, 'Position', [00, 00, 400, 300])
set(gca,'xtick',[0,125,250,375,500]);
set(gca,'xticklabel',{'-0.5','-0.25','0','0.25','0.5'});
axis([-0 N 0 1.1])

figure
plot(real(Ksi),'b*')
hold on
plot(imag(Ksi),'r*')
set(gca,'FontSize', 14)
set(gcf, 'Position', [00, 00, 400, 300])
set(gca,'xtick',[0,N/4,N/2,N/4*3,N]);
set(gca,'xticklabel',{'-0.5','-0.25','0','0.25','0.5'});
axis([0 N -2 12])

figure
plot(real(n),'b*')
hold on
plot(imag(n),'r*')
set(gca,'FontSize', 14)
set(gcf, 'Position', [00, 00, 400, 300])
set(gca,'xtick',[0,N/4,N/2,N/4*3,N])
set(gca,'xticklabel',{'-0.5','-0.25','0','0.25','0.5'});
axis([0 N -1 3.5])

yt = sign(x);
idx = yt==0;
yt(idx) = 1;
wx = yt.*real(n);
figure
plot(Sx,'b-')
set(gca,'FontSize', 14)
set(gcf, 'Position', [00, 00, 400, 300])
set(gca,'xtick',[0,N/4,N/2,N/4*3,N])
set(gca,'xticklabel',{'-0.5','-0.25','0','0.25','0.5'});
axis([0 N -3.5 3.5])