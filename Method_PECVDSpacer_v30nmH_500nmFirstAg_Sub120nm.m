clc
clear all
close all

Data = [1.51000000000000e-07 + 0.00000000000000i,0.0151229830536482 + 0.189522549513153i;1.52000000000000e-07 + 0.00000000000000i,0.209511737072606 + 0.0155512821299974i;1.53000000000000e-07 + 0.00000000000000i,0.349402731746243 + 0.0104134089163309i;1.54000000000000e-07 + 0.00000000000000i,0.446492362395337 + 0.00897500235648818i;1.55000000000000e-07 + 0.00000000000000i,0.525001382354219 + 0.00831256979566498i;1.56000000000000e-07 + 0.00000000000000i,0.592398492229176 + 0.00795097082492759i;1.57000000000000e-07 + 0.00000000000000i,0.652141995370397 + 0.00773707380057277i;1.58000000000000e-07 + 0.00000000000000i,0.706131216089037 + 0.00760262065434706i;1.59000000000000e-07 + 0.00000000000000i,0.755677438174187 + 0.00751925699505842i;1.60000000000000e-07 + 0.00000000000000i,0.801580266047164 + 0.00746672063762601i;1.61000000000000e-07 + 0.00000000000000i,0.844407568808611 + 0.00743633536239645i;1.62000000000000e-07 + 0.00000000000000i,0.884617108971818 + 0.00741997172522361i;1.63000000000000e-07 + 0.00000000000000i,0.922554750595572 + 0.00740839303040287i;1.64000000000000e-07 + 0.00000000000000i,0.958438544941201 + 0.00738417050280933i;1.65000000000000e-07 + 0.00000000000000i,0.992443256736036 + 0.00733071424480431i;1.66000000000000e-07 + 0.00000000000000i,1.02469029492002 + 0.00722425181719741i;1.67000000000000e-07 + 0.00000000000000i,1.05519854887894 + 0.00701865029186107i;1.68000000000000e-07 + 0.00000000000000i,1.08383132157697 + 0.00663503852790920i;1.69000000000000e-07 + 0.00000000000000i,1.11033111195978 + 0.00600472882842459i;1.70000000000000e-07 + 0.00000000000000i,1.15635079270543 + 0.00587248495872299i;1.71000000000000e-07 + 0.00000000000000i,1.18088446977690 + 0.00650800091753958i;1.72000000000000e-07 + 0.00000000000000i,1.20549949755267 + 0.00692660468160499i;1.73000000000000e-07 + 0.00000000000000i,1.22983672318909 + 0.00718912107321660i;1.74000000000000e-07 + 0.00000000000000i,1.25370911601775 + 0.00736062240038960i;1.75000000000000e-07 + 0.00000000000000i,1.27702330443997 + 0.00748002363072209i;1.76000000000000e-07 + 0.00000000000000i,1.29976669618545 + 0.00756912056403154i;1.77000000000000e-07 + 0.00000000000000i,1.32192354443051 + 0.00763926024239887i;1.78000000000000e-07 + 0.00000000000000i,1.34352304352141 + 0.00769730659313413i;1.79000000000000e-07 + 0.00000000000000i,1.36457300393093 + 0.00774709491611727i;1.80000000000000e-07 + 0.00000000000000i,1.38508460704893 + 0.00779078087308236i;1.81000000000000e-07 + 0.00000000000000i,1.40510456636823 + 0.00783038694759104i;1.82000000000000e-07 + 0.00000000000000i,1.42463664434925 + 0.00786690100954007i;1.83000000000000e-07 + 0.00000000000000i,1.44370354825151 + 0.00790049534849058i;1.84000000000000e-07 + 0.00000000000000i,1.46231716224374 + 0.00793172438911797i;1.85000000000000e-07 + 0.00000000000000i,1.48051443482192 + 0.00796154247569220i;1.86000000000000e-07 + 0.00000000000000i,1.49829206475004 + 0.00798962684507071i;1.87000000000000e-07 + 0.00000000000000i,1.51568249010809 + 0.00801642200202350i;1.88000000000000e-07 + 0.00000000000000i,1.53268984032168 + 0.00804221708123226i;1.89000000000000e-07 + 0.00000000000000i,1.54933542552535 + 0.00806679039337684i;1.90000000000000e-07 + 0.00000000000000i,1.56563387627225 + 0.00809055575288051i;1.91000000000000e-07 + 0.00000000000000i,1.58158962558325 + 0.00811332532458612i;1.92000000000000e-07 + 0.00000000000000i,1.59722086536689 + 0.00813530721129641i;1.93000000000000e-07 + 0.00000000000000i,1.61254288642829 + 0.00815671997858444i;1.94000000000000e-07 + 0.00000000000000i,1.62755952176536 + 0.00817728202833598i;1.95000000000000e-07 + 0.00000000000000i,1.64229161760085 + 0.00819712045785445i;1.96000000000000e-07 + 0.00000000000000i,1.65673600807835 + 0.00821643596547288i;1.97000000000000e-07 + 0.00000000000000i,1.67090692365814 + 0.00823525982610220i;1.98000000000000e-07 + 0.00000000000000i,1.68481330973845 + 0.00825333787137608i;1.99000000000000e-07 + 0.00000000000000i,1.69845958561985 + 0.00827092090968147i;2.00000000000000e-07 + 0.00000000000000i,1.71187280596925 + 0.00828831902238274i;2.01000000000000e-07 + 0.00000000000000i,1.72503971420583 + 0.00830505005192264i;2.02000000000000e-07 + 0.00000000000000i,1.73796659093919 + 0.00832111222560393i;2.03000000000000e-07 + 0.00000000000000i,1.75067973215420 + 0.00833712727910752i;2.04000000000000e-07 + 0.00000000000000i,1.76316441847604 + 0.00835239674568422i;2.05000000000000e-07 + 0.00000000000000i,1.77543949909990 + 0.00836735409492009i;2.06000000000000e-07 + 0.00000000000000i,1.78750598157428 + 0.00838193935939396i;2.07000000000000e-07 + 0.00000000000000i,1.79937718091947 + 0.00839621919537095i;2.08000000000000e-07 + 0.00000000000000i,1.81104814221781 + 0.00840999867890673i;2.09000000000000e-07 + 0.00000000000000i,1.82253194193165 + 0.00842349421894338i;2.10000000000000e-07 + 0.00000000000000i,1.83382771124927 + 0.00843652949828351i;2.11000000000000e-07 + 0.00000000000000i,1.84494755191896 + 0.00844937046465363i;2.12000000000000e-07 + 0.00000000000000i,1.85588722242739 + 0.00846169004489125i;2.13000000000000e-07 + 0.00000000000000i,1.86666133224008 + 0.00847391145669590i;2.14000000000000e-07 + 0.00000000000000i,1.87726318762294 + 0.00848565550637034i;2.15000000000000e-07 + 0.00000000000000i,1.88771061449125 + 0.00849732604282147i;2.16000000000000e-07 + 0.00000000000000i,1.89798768183101 + 0.00850840583063182i;2.17000000000000e-07 + 0.00000000000000i,1.90811882157316 + 0.00851941053914205i;2.18000000000000e-07 + 0.00000000000000i,1.91809997908088 + 0.00853015279591855i;2.19000000000000e-07 + 0.00000000000000i,1.92793037684870 + 0.00854049932436848i;2.20000000000000e-07 + 0.00000000000000i,1.93761710168734 + 0.00855068463375317i;2.21000000000000e-07 + 0.00000000000000i,1.94717244972654 + 0.00856069117703894i;2.22000000000000e-07 + 0.00000000000000i,1.95658394767659 + 0.00857034847056768i;2.23000000000000e-07 + 0.00000000000000i,1.96586311463217 + 0.00857986742046894i;2.24000000000000e-07 + 0.00000000000000i,1.97500868591597 + 0.00858914495750659i;2.25000000000000e-07 + 0.00000000000000i,1.98402650272658 + 0.00859811345580399i;2.26000000000000e-07 + 0.00000000000000i,1.99292199961338 + 0.00860700821088722i;2.27000000000000e-07 + 0.00000000000000i,2.00168766924803 + 0.00861540497813530i;2.28000000000000e-07 + 0.00000000000000i,2.01033888973673 + 0.00862380010569043i;2.29000000000000e-07 + 0.00000000000000i,2.01887217964812 + 0.00863202226012528i;2.30000000000000e-07 + 0.00000000000000i,2.02729139499598 + 0.00864006535040346i;2.32000000000000e-07 + 0.00000000000000i,2.04379385205901 + 0.00865555855630893i;2.34000000000000e-07 + 0.00000000000000i,2.05986578469701 + 0.00867021772627012i;2.36000000000000e-07 + 0.00000000000000i,2.07552067141405 + 0.00868428196220798i;2.38000000000000e-07 + 0.00000000000000i,2.09078088211021 + 0.00869765926143533i;2.40000000000000e-07 + 0.00000000000000i,2.10565587139292 + 0.00871042933466028i;2.42000000000000e-07 + 0.00000000000000i,2.12016319709754 + 0.00872252068671461i;2.44000000000000e-07 + 0.00000000000000i,2.13431972136181 + 0.00873433693249614i;2.46000000000000e-07 + 0.00000000000000i,2.14813289937546 + 0.00874524529952080i;2.48000000000000e-07 + 0.00000000000000i,2.16161895199830 + 0.00875591253380863i;2.50000000000000e-07 + 0.00000000000000i,2.17478551480569 + 0.00876603886417263i;2.52000000000000e-07 + 0.00000000000000i,2.18764787478993 + 0.00877573012325666i;2.54000000000000e-07 + 0.00000000000000i,2.20021227361782 + 0.00878496548472628i;2.56000000000000e-07 + 0.00000000000000i,2.21248962832870 + 0.00879371340759106i;2.58000000000000e-07 + 0.00000000000000i,2.22449292056747 + 0.00880208284108151i;2.60000000000000e-07 + 0.00000000000000i,2.23623235168660 + 0.00881017901895016i;2.62000000000000e-07 + 0.00000000000000i,2.24771219910080 + 0.00881782772035769i;2.64000000000000e-07 + 0.00000000000000i,2.25894179109751 + 0.00882509508460186i;2.66000000000000e-07 + 0.00000000000000i,2.26992996846068 + 0.00883205108964773i;2.68000000000000e-07 + 0.00000000000000i,2.28068224220484 + 0.00883867349061183i;2.70000000000000e-07 + 0.00000000000000i,2.29120922607967 + 0.00884499974408502i;2.72000000000000e-07 + 0.00000000000000i,2.30140117791357 + 0.00884769496612726i;2.74000000000000e-07 + 0.00000000000000i,2.31150046821777 + 0.00885363142109189i;2.76000000000000e-07 + 0.00000000000000i,2.32139677097587 + 0.00885928174817029i;2.78000000000000e-07 + 0.00000000000000i,2.33108267162145 + 0.00886452053357876i;2.80000000000000e-07 + 0.00000000000000i,2.34057721952478 + 0.00886957037017093i;2.82000000000000e-07 + 0.00000000000000i,2.34988250085819 + 0.00887435617394465i;2.84000000000000e-07 + 0.00000000000000i,2.35899751392551 + 0.00887876146816037i;2.86000000000000e-07 + 0.00000000000000i,2.36794415514656 + 0.00888323087086264i;2.88000000000000e-07 + 0.00000000000000i,2.37671349714321 + 0.00888735142051091i;2.90000000000000e-07 + 0.00000000000000i,2.38531226778390 + 0.00889123279273514i;2.92000000000000e-07 + 0.00000000000000i,2.39374860222649 + 0.00889499463156288i;2.94000000000000e-07 + 0.00000000000000i,2.40202611707127 + 0.00889855192815246i;2.96000000000000e-07 + 0.00000000000000i,2.41014433409039 + 0.00890180501693385i;2.98000000000000e-07 + 0.00000000000000i,2.41811534965920 + 0.00890501412430239i;3.00000000000000e-07 + 0.00000000000000i,2.42593640456263 + 0.00890797383810105i;3.02000000000000e-07 + 0.00000000000000i,2.43361690702525 + 0.00891095257070115i;3.04000000000000e-07 + 0.00000000000000i,2.44116057235968 + 0.00891369535393456i;3.06000000000000e-07 + 0.00000000000000i,2.44855985836099 + 0.00891610980247352i;3.08000000000000e-07 + 0.00000000000000i,2.45583498991445 + 0.00891856524884417i;3.10000000000000e-07 + 0.00000000000000i,2.46297879579139 + 0.00892090935870884i;3.12000000000000e-07 + 0.00000000000000i,2.46999505567218 + 0.00892304729449396i;3.14000000000000e-07 + 0.00000000000000i,2.47688182113878 + 0.00892483900832744i;3.16000000000000e-07 + 0.00000000000000i,2.48366021894473 + 0.00892686154113294i;3.18000000000000e-07 + 0.00000000000000i,2.49032152612865 + 0.00892878420352856i;3.20000000000000e-07 + 0.00000000000000i,2.49686140101486 + 0.00893033875448090i;3.22000000000000e-07 + 0.00000000000000i,2.50329101886185 + 0.00893183225352998i;3.24000000000000e-07 + 0.00000000000000i,2.50961718412090 + 0.00893335413063876i;3.26000000000000e-07 + 0.00000000000000i,2.51583389745970 + 0.00893471293575899i;3.28000000000000e-07 + 0.00000000000000i,2.52194990726258 + 0.00893603622726250i;3.30000000000000e-07 + 0.00000000000000i,2.52795874643605 + 0.00893707811350933i;3.32000000000000e-07 + 0.00000000000000i,2.53387405542325 + 0.00893821043079961i;3.34000000000000e-07 + 0.00000000000000i,2.53969178089240 + 0.00893919707649587i;3.36000000000000e-07 + 0.00000000000000i,2.54541648835690 + 0.00894007238666919i;3.38000000000000e-07 + 0.00000000000000i,2.55104503712652 + 0.00894085311336756i;3.40000000000000e-07 + 0.00000000000000i,2.55658778441060 + 0.00894164484597327i;3.42000000000000e-07 + 0.00000000000000i,2.56204087109698 + 0.00894230127719996i;3.44000000000000e-07 + 0.00000000000000i,2.56740685472507 + 0.00894293858187132i;3.46000000000000e-07 + 0.00000000000000i,2.57268599665028 + 0.00894336787051752i;3.48000000000000e-07 + 0.00000000000000i,2.57788989072575 + 0.00894394114940490i;3.50000000000000e-07 + 0.00000000000000i,2.58300957781256 + 0.00894434120556401i;3.52000000000000e-07 + 0.00000000000000i,2.58804885346789 + 0.00894467916692421i;3.54000000000000e-07 + 0.00000000000000i,2.59300807746808 + 0.00894480850961809i;3.56000000000000e-07 + 0.00000000000000i,2.59790068795356 + 0.00894520829212456i;3.58000000000000e-07 + 0.00000000000000i,2.60271433752733 + 0.00894536236223014i;3.60000000000000e-07 + 0.00000000000000i,2.60745538589967 + 0.00894546504296555i;3.62000000000000e-07 + 0.00000000000000i,2.61212539335484 + 0.00894550267819151i;3.64000000000000e-07 + 0.00000000000000i,2.61672713910652 + 0.00894555429602519i;3.66000000000000e-07 + 0.00000000000000i,2.62125831464775 + 0.00894547428353826i;3.68000000000000e-07 + 0.00000000000000i,2.62572882581788 + 0.00894547337334459i;3.70000000000000e-07 + 0.00000000000000i,2.63012910909995 + 0.00894533095756269i;3.72000000000000e-07 + 0.00000000000000i,2.63446602225623 + 0.00894517488631244i;3.74000000000000e-07 + 0.00000000000000i,2.63874058762943 + 0.00894495895255838i;3.76000000000000e-07 + 0.00000000000000i,2.64295200620612 + 0.00894469419539013i;3.78000000000000e-07 + 0.00000000000000i,2.64710663399542 + 0.00894446445120216i;3.80000000000000e-07 + 0.00000000000000i,2.65119757511723 + 0.00894407849444111i;3.82000000000000e-07 + 0.00000000000000i,2.65523381056886 + 0.00894375066873915i;3.84000000000000e-07 + 0.00000000000000i,2.65921441044036 + 0.00894341690876518i;3.86000000000000e-07 + 0.00000000000000i,2.66313704204461 + 0.00894298026081595i;3.88000000000000e-07 + 0.00000000000000i,2.66700546544549 + 0.00894252513851714i;3.90000000000000e-07 + 0.00000000000000i,2.67081948769982 + 0.00894199011903750i;3.92000000000000e-07 + 0.00000000000000i,2.67458477115687 + 0.00894156856192335i;3.94000000000000e-07 + 0.00000000000000i,2.67829435953409 + 0.00894095416339519i;3.96000000000000e-07 + 0.00000000000000i,2.68195464273734 + 0.00894040006130262i;3.98000000000000e-07 + 0.00000000000000i,2.68556699476122 + 0.00893982959096968i;4.00000000000000e-07 + 0.00000000000000i,2.68913018140267 + 0.00893923117814997i];
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
title('120nm InP Left');
xlim([130 500]);

k0 = 2*pi/1550e-5;
L = 20e-6;
M = 300;
ne = 2.6;
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
SiO_Y = [d_Al + 60e-9 -d_Al - 60e-9]/2*1e6;
SiO_Z = [z -z];


% % %%%%3.3^2 +(0.03*exp(-((z-5.168[um])/2[um])^2) + 0.03*exp(-((z+5.168[um])/2[um])^2) + 0.017*exp(-(z/10[um])^2))*i
% % tm = 0.03*exp(-((z-5.168e-6)/2e-6).^2) + 0.03*exp(-((z+5.168e-6)/2e-6).^2) + 0.017*exp(-(z/10e-6).^2);
% % figure
% % plot(z,tm)