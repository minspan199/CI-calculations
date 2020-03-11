clc
clear all
close all

Data = [1.47000000000000e-07 + 0.00000000000000i,0.00459265375349053 + 0.279833593934654i;1.48000000000000e-07 + 0.00000000000000i,0.0736334213031090 + 0.0194157008755764i;1.49000000000000e-07 + 0.00000000000000i,0.295332817561764 + 0.00557107762332136i;1.50000000000000e-07 + 0.00000000000000i,0.410088790560833 + 0.00459864115284085i;1.51000000000000e-07 + 0.00000000000000i,0.497959983786653 + 0.00425820916523453i;1.52000000000000e-07 + 0.00000000000000i,0.571490946517607 + 0.00406819474066138i;1.53000000000000e-07 + 0.00000000000000i,0.635742322209391 + 0.00389626664073760i;1.54000000000000e-07 + 0.00000000000000i,0.693355028201377 + 0.00367629672271683i;1.55000000000000e-07 + 0.00000000000000i,0.745955941424516 + 0.00336488391715683i;1.56000000000000e-07 + 0.00000000000000i,0.794766360697258 + 0.00299628996350575i;1.57000000000000e-07 + 0.00000000000000i,0.840756718560887 + 0.00303775611313557i;1.58000000000000e-07 + 0.00000000000000i,0.883457187900586 + 0.00425450871840439i;1.59000000000000e-07 + 0.00000000000000i,0.922394362204785 + 0.00531393501765413i;1.60000000000000e-07 + 0.00000000000000i,0.958713382538675 + 0.00558176527091828i;1.61000000000000e-07 + 0.00000000000000i,0.993018906897684 + 0.00531589363444733i;1.62000000000000e-07 + 0.00000000000000i,1.02526052838445 + 0.00458129374353243i;1.63000000000000e-07 + 0.00000000000000i,1.05521383163329 + 0.00327461393093802i;1.64000000000000e-07 + 0.00000000000000i,1.08234386019758 + 0.00151150700726028i;1.65000000000000e-07 + 0.00000000000000i,1.13569587363761 + 0.000743301610814812i;1.66000000000000e-07 + 0.00000000000000i,1.16253494754153 + 0.00269368008258892i;1.67000000000000e-07 + 0.00000000000000i,1.18910262546402 + 0.00386128476695415i;1.68000000000000e-07 + 0.00000000000000i,1.21512838508889 + 0.00458608228868267i;1.69000000000000e-07 + 0.00000000000000i,1.24051672802453 + 0.00506222437278725i;1.70000000000000e-07 + 0.00000000000000i,1.26525096978817 + 0.00539324286042997i;1.71000000000000e-07 + 0.00000000000000i,1.28931107650137 + 0.00563494606208103i;1.72000000000000e-07 + 0.00000000000000i,1.31273330408996 + 0.00582002358240575i;1.73000000000000e-07 + 0.00000000000000i,1.33551534073020 + 0.00596650333240594i;1.74000000000000e-07 + 0.00000000000000i,1.35770442603680 + 0.00608697877058464i;1.75000000000000e-07 + 0.00000000000000i,1.37931139423960 + 0.00618824115724777i;1.76000000000000e-07 + 0.00000000000000i,1.40036743089210 + 0.00627550451595759i;1.77000000000000e-07 + 0.00000000000000i,1.42088798490206 + 0.00635196498878481i;1.78000000000000e-07 + 0.00000000000000i,1.44090902505673 + 0.00642023111333383i;1.79000000000000e-07 + 0.00000000000000i,1.46045007096497 + 0.00648206306006900i;1.80000000000000e-07 + 0.00000000000000i,1.47950923969852 + 0.00653805385155160i;1.81000000000000e-07 + 0.00000000000000i,1.49814054573231 + 0.00659004497154295i;1.82000000000000e-07 + 0.00000000000000i,1.51634752847710 + 0.00663853048517401i;1.83000000000000e-07 + 0.00000000000000i,1.53414044551157 + 0.00668339887476261i;1.84000000000000e-07 + 0.00000000000000i,1.55154339566582 + 0.00672578850242535i;1.85000000000000e-07 + 0.00000000000000i,1.56857526448151 + 0.00676594507815394i;1.86000000000000e-07 + 0.00000000000000i,1.58523761970470 + 0.00680379641573166i;1.87000000000000e-07 + 0.00000000000000i,1.60155667683545 + 0.00683985737693934i;1.88000000000000e-07 + 0.00000000000000i,1.61753779514687 + 0.00687435972694960i;1.89000000000000e-07 + 0.00000000000000i,1.63319620113113 + 0.00690723001792440i;1.90000000000000e-07 + 0.00000000000000i,1.64854604748194 + 0.00693886272670667i;1.91000000000000e-07 + 0.00000000000000i,1.66358619695067 + 0.00696901083742458i;1.92000000000000e-07 + 0.00000000000000i,1.67834152545528 + 0.00699817315838323i;1.93000000000000e-07 + 0.00000000000000i,1.69281693021054 + 0.00702633947105031i;1.94000000000000e-07 + 0.00000000000000i,1.70701688763881 + 0.00705337690212012i;1.95000000000000e-07 + 0.00000000000000i,1.72095968473161 + 0.00707970408861855i;1.96000000000000e-07 + 0.00000000000000i,1.73464254928064 + 0.00710473436477296i;1.97000000000000e-07 + 0.00000000000000i,1.74807789934797 + 0.00712922123049894i;1.98000000000000e-07 + 0.00000000000000i,1.76127374150084 + 0.00715281045685773i;1.99000000000000e-07 + 0.00000000000000i,1.77423432893371 + 0.00717565235516055i;2.00000000000000e-07 + 0.00000000000000i,1.78697992973025 + 0.00719812339703034i;2.01000000000000e-07 + 0.00000000000000i,1.79950173642891 + 0.00721976355776602i;2.02000000000000e-07 + 0.00000000000000i,1.81180770860263 + 0.00724067181848609i;2.03000000000000e-07 + 0.00000000000000i,1.82391049133088 + 0.00726113086352537i;2.04000000000000e-07 + 0.00000000000000i,1.83580932912038 + 0.00728087743928186i;2.05000000000000e-07 + 0.00000000000000i,1.84751436603970 + 0.00730014236781366i;2.06000000000000e-07 + 0.00000000000000i,1.85902783291659 + 0.00731888374251212i;2.07000000000000e-07 + 0.00000000000000i,1.87036145735601 + 0.00733722333833724i;2.08000000000000e-07 + 0.00000000000000i,1.88150984853993 + 0.00735494410694677i;2.09000000000000e-07 + 0.00000000000000i,1.89248618188831 + 0.00737228403805156i;2.10000000000000e-07 + 0.00000000000000i,1.90328991605977 + 0.00738912354132534i;2.11000000000000e-07 + 0.00000000000000i,1.91392519651145 + 0.00740553489625320i;2.12000000000000e-07 + 0.00000000000000i,1.92439732538720 + 0.00742143866473667i;2.13000000000000e-07 + 0.00000000000000i,1.93471506388374 + 0.00743712043728049i;2.14000000000000e-07 + 0.00000000000000i,1.94487335414890 + 0.00745224975655315i;2.15000000000000e-07 + 0.00000000000000i,1.95488765334342 + 0.00746723299260955i;2.16000000000000e-07 + 0.00000000000000i,1.96474497232188 + 0.00748153570281153i;2.17000000000000e-07 + 0.00000000000000i,1.97446418382471 + 0.00749570334089257i;2.18000000000000e-07 + 0.00000000000000i,1.98404373621525 + 0.00750954499737814i;2.19000000000000e-07 + 0.00000000000000i,1.99348207296975 + 0.00752295672389890i;2.20000000000000e-07 + 0.00000000000000i,2.00278687077827 + 0.00753606139211174i;2.21000000000000e-07 + 0.00000000000000i,2.01196790132708 + 0.00754901177598872i;2.22000000000000e-07 + 0.00000000000000i,2.02101289216103 + 0.00756152272422216i;2.23000000000000e-07 + 0.00000000000000i,2.02993730322918 + 0.00757384815941059i;2.24000000000000e-07 + 0.00000000000000i,2.03873471593211 + 0.00758584566755618i;2.25000000000000e-07 + 0.00000000000000i,2.04741249018366 + 0.00759752486229406i;2.26000000000000e-07 + 0.00000000000000i,2.05597496934214 + 0.00760904525767278i;2.27000000000000e-07 + 0.00000000000000i,2.06441576786503 + 0.00762004410248266i;2.28000000000000e-07 + 0.00000000000000i,2.07274799539582 + 0.00763098349715005i;2.29000000000000e-07 + 0.00000000000000i,2.08097071859198 + 0.00764170478066250i;2.30000000000000e-07 + 0.00000000000000i,2.08908460140491 + 0.00765219249146529i;2.32000000000000e-07 + 0.00000000000000i,2.10499689895027 + 0.00767246447605363i;2.34000000000000e-07 + 0.00000000000000i,2.12050149570448 + 0.00769177804257948i;2.36000000000000e-07 + 0.00000000000000i,2.13561263422592 + 0.00771030462776957i;2.38000000000000e-07 + 0.00000000000000i,2.15035029739825 + 0.00772803731300384i;2.40000000000000e-07 + 0.00000000000000i,2.16472303306827 + 0.00774500920761633i;2.42000000000000e-07 + 0.00000000000000i,2.17874650399898 + 0.00776120384846228i;2.44000000000000e-07 + 0.00000000000000i,2.19243729907239 + 0.00777695128070626i;2.46000000000000e-07 + 0.00000000000000i,2.20580062093223 + 0.00779177044269776i;2.48000000000000e-07 + 0.00000000000000i,2.21885396699452 + 0.00780617335961622i;2.50000000000000e-07 + 0.00000000000000i,2.23160250865630 + 0.00781994463697695i;2.52000000000000e-07 + 0.00000000000000i,2.24406178920780 + 0.00783318207965615i;2.54000000000000e-07 + 0.00000000000000i,2.25623730612035 + 0.00784583534184683i;2.56000000000000e-07 + 0.00000000000000i,2.26813902830287 + 0.00785794154589326i;2.58000000000000e-07 + 0.00000000000000i,2.27977835952501 + 0.00786958015489868i;2.60000000000000e-07 + 0.00000000000000i,2.29116538875257 + 0.00788085332865972i;2.62000000000000e-07 + 0.00000000000000i,2.30230395660881 + 0.00789160886478961i;2.64000000000000e-07 + 0.00000000000000i,2.31320341932169 + 0.00790191549654715i;2.66000000000000e-07 + 0.00000000000000i,2.32387145584526 + 0.00791182751056502i;2.68000000000000e-07 + 0.00000000000000i,2.33431372566274 + 0.00792132763465695i;2.70000000000000e-07 + 0.00000000000000i,2.34453964574651 + 0.00793046923172784i;2.72000000000000e-07 + 0.00000000000000i,2.35455818306806 + 0.00793930799333832i;2.74000000000000e-07 + 0.00000000000000i,2.36436974205391 + 0.00794776861070118i;2.76000000000000e-07 + 0.00000000000000i,2.37397927090522 + 0.00795579442939875i;2.78000000000000e-07 + 0.00000000000000i,2.38339924420913 + 0.00796360748032480i;2.80000000000000e-07 + 0.00000000000000i,2.39263098915335 + 0.00797106883053444i;2.82000000000000e-07 + 0.00000000000000i,2.40168049637503 + 0.00797823734549423i;2.84000000000000e-07 + 0.00000000000000i,2.41055479751367 + 0.00798520165113345i;2.86000000000000e-07 + 0.00000000000000i,2.41925721088651 + 0.00799187466350456i;2.88000000000000e-07 + 0.00000000000000i,2.42778964958585 + 0.00799818965790504i;2.90000000000000e-07 + 0.00000000000000i,2.43616279566879 + 0.00800433499143316i;2.92000000000000e-07 + 0.00000000000000i,2.44437530450201 + 0.00801017486779941i;2.94000000000000e-07 + 0.00000000000000i,2.45243473557513 + 0.00801581491472929i;2.96000000000000e-07 + 0.00000000000000i,2.46034574775150 + 0.00802122970470426i;2.98000000000000e-07 + 0.00000000000000i,2.46810956802437 + 0.00802643945472703i;3.00000000000000e-07 + 0.00000000000000i,2.47573272723644 + 0.00803145959251565i;3.02000000000000e-07 + 0.00000000000000i,2.48321638777015 + 0.00803625681282123i;3.04000000000000e-07 + 0.00000000000000i,2.49056627073355 + 0.00804089434750967i;3.06000000000000e-07 + 0.00000000000000i,2.49778402272479 + 0.00804531992683752i;3.08000000000000e-07 + 0.00000000000000i,2.50487393266386 + 0.00804957496500394i;3.10000000000000e-07 + 0.00000000000000i,2.51183877034371 + 0.00805366325983253i;3.12000000000000e-07 + 0.00000000000000i,2.51868187272003 + 0.00805756745175661i;3.14000000000000e-07 + 0.00000000000000i,2.52540565294325 + 0.00806131429977955i;3.16000000000000e-07 + 0.00000000000000i,2.53201545380448 + 0.00806493341216811i;3.18000000000000e-07 + 0.00000000000000i,2.53851036363188 + 0.00806835270479844i;3.20000000000000e-07 + 0.00000000000000i,2.54489729635137 + 0.00807169101730553i;3.22000000000000e-07 + 0.00000000000000i,2.55117374202497 + 0.00807483045556440i;3.24000000000000e-07 + 0.00000000000000i,2.55734611882951 + 0.00807786068731850i;3.26000000000000e-07 + 0.00000000000000i,2.56341633576098 + 0.00808078279471464i;3.28000000000000e-07 + 0.00000000000000i,2.56938574360220 + 0.00808356020225960i;3.30000000000000e-07 + 0.00000000000000i,2.57525609528402 + 0.00808620660369399i;3.32000000000000e-07 + 0.00000000000000i,2.58103195758686 + 0.00808875816702861i;3.34000000000000e-07 + 0.00000000000000i,2.58671446602341 + 0.00809120844784927i;3.36000000000000e-07 + 0.00000000000000i,2.59230458897460 + 0.00809351792616044i;3.38000000000000e-07 + 0.00000000000000i,2.59780588911748 + 0.00809575256248255i;3.40000000000000e-07 + 0.00000000000000i,2.60321990872517 + 0.00809788729846111i;3.42000000000000e-07 + 0.00000000000000i,2.60854787128959 + 0.00809991168789086i;3.44000000000000e-07 + 0.00000000000000i,2.61379273205583 + 0.00810184483381109i;3.46000000000000e-07 + 0.00000000000000i,2.61895612163663 + 0.00810369440430487i;3.48000000000000e-07 + 0.00000000000000i,2.62403841123180 + 0.00810542132287442i;3.50000000000000e-07 + 0.00000000000000i,2.62904326658563 + 0.00810711071776788i;3.52000000000000e-07 + 0.00000000000000i,2.63397033642497 + 0.00810867356037281i;3.54000000000000e-07 + 0.00000000000000i,2.63882405262962 + 0.00811021142349981i;3.56000000000000e-07 + 0.00000000000000i,2.64360335307075 + 0.00811164010333704i;3.58000000000000e-07 + 0.00000000000000i,2.64831006495338 + 0.00811298758101153i;3.60000000000000e-07 + 0.00000000000000i,2.65294762720608 + 0.00811428147371353i;3.62000000000000e-07 + 0.00000000000000i,2.65751560263535 + 0.00811549828866977i;3.64000000000000e-07 + 0.00000000000000i,2.66201558241647 + 0.00811662130623811i;3.66000000000000e-07 + 0.00000000000000i,2.66645058352608 + 0.00811773062542897i;3.68000000000000e-07 + 0.00000000000000i,2.67082062905421 + 0.00811876322337173i;3.70000000000000e-07 + 0.00000000000000i,2.67512612165450 + 0.00811969327591659i;3.72000000000000e-07 + 0.00000000000000i,2.67937069400796 + 0.00812061562177164i;3.74000000000000e-07 + 0.00000000000000i,2.68355351534732 + 0.00812146033773958i;3.76000000000000e-07 + 0.00000000000000i,2.68767627757135 + 0.00812225412351786i;3.78000000000000e-07 + 0.00000000000000i,2.69174066035094 + 0.00812299301569475i;3.80000000000000e-07 + 0.00000000000000i,2.69574771282408 + 0.00812367742659718i;3.82000000000000e-07 + 0.00000000000000i,2.69969825822914 + 0.00812431577039215i;3.84000000000000e-07 + 0.00000000000000i,2.70359350592985 + 0.00812490986951258i;3.86000000000000e-07 + 0.00000000000000i,2.70743503636954 + 0.00812546001352128i;3.88000000000000e-07 + 0.00000000000000i,2.71122293950673 + 0.00812596028985001i;3.90000000000000e-07 + 0.00000000000000i,2.71495856625593 + 0.00812640916356151i;3.92000000000000e-07 + 0.00000000000000i,2.71864372598377 + 0.00812684350596592i;3.94000000000000e-07 + 0.00000000000000i,2.72227805928393 + 0.00812720899199267i;3.96000000000000e-07 + 0.00000000000000i,2.72586316264289 + 0.00812752927808316i;3.98000000000000e-07 + 0.00000000000000i,2.72940097943588 + 0.00812785728505292i;4.00000000000000e-07 + 0.00000000000000i,2.73288950567394 + 0.00812809287364256i];
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
title('250nm InP Left');
xlim([130 500]);
% legend(["50nm","100nm","120nm","150nm","200nm","250nm"])

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