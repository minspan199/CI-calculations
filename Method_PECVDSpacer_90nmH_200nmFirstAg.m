clc
clear all
% close all

% Htop = 3*Ht
Data = [1.07000000000000e-07 + 0.00000000000000i,0.0271053774187246 + 0.168296866325157i;1.08000000000000e-07 + 0.00000000000000i,0.0889074336078656 + 0.0508667000882510i;1.09000000000000e-07 + 0.00000000000000i,0.195684112126023 + 0.0229113942203471i;1.10000000000000e-07 + 0.00000000000000i,0.264674043262599 + 0.0167929419115399i;1.11000000000000e-07 + 0.00000000000000i,0.318646198135128 + 0.0138281026709676i;1.12000000000000e-07 + 0.00000000000000i,0.364183236528540 + 0.0119941066184203i;1.13000000000000e-07 + 0.00000000000000i,0.404128062850030 + 0.0107142865060655i;1.14000000000000e-07 + 0.00000000000000i,0.440021950952493 + 0.00975444351616230i;1.15000000000000e-07 + 0.00000000000000i,0.472780194021541 + 0.00899885966991101i;1.16000000000000e-07 + 0.00000000000000i,0.503018736289634 + 0.00838296222509983i;1.17000000000000e-07 + 0.00000000000000i,0.531184565988212 + 0.00786770466486579i;1.18000000000000e-07 + 0.00000000000000i,0.557595830430897 + 0.00742768164821658i;1.19000000000000e-07 + 0.00000000000000i,0.582492835108862 + 0.00704551029749124i;1.20000000000000e-07 + 0.00000000000000i,0.606072117229495 + 0.00670908957253196i;1.21000000000000e-07 + 0.00000000000000i,0.628488646663604 + 0.00640946494943275i;1.22000000000000e-07 + 0.00000000000000i,0.649867348772651 + 0.00613990441389302i;1.23000000000000e-07 + 0.00000000000000i,0.670308469393215 + 0.00589516122526598i;1.24000000000000e-07 + 0.00000000000000i,0.689899574593070 + 0.00567113060091464i;1.25000000000000e-07 + 0.00000000000000i,0.708725309867822 + 0.00546485866382507i;1.26000000000000e-07 + 0.00000000000000i,0.726840592857655 + 0.00527348162381885i;1.27000000000000e-07 + 0.00000000000000i,0.744303185996220 + 0.00509487609170853i;1.28000000000000e-07 + 0.00000000000000i,0.761167497409794 + 0.00492735962259548i;1.29000000000000e-07 + 0.00000000000000i,0.777471680640955 + 0.00476933287732910i;1.30000000000000e-07 + 0.00000000000000i,0.793255581562548 + 0.00461955645412963i;1.31000000000000e-07 + 0.00000000000000i,0.808552809139217 + 0.00447698476711290i;1.32000000000000e-07 + 0.00000000000000i,0.823390998339826 + 0.00434064669385193i;1.33000000000000e-07 + 0.00000000000000i,0.837804480398429 + 0.00421003236290581i;1.34000000000000e-07 + 0.00000000000000i,0.851814019961251 + 0.00408452177778057i;1.35000000000000e-07 + 0.00000000000000i,0.865444734513892 + 0.00396380785644012i;1.36000000000000e-07 + 0.00000000000000i,0.878717557028859 + 0.00384779557596921i;1.37000000000000e-07 + 0.00000000000000i,0.891651329554199 + 0.00373652684285329i;1.38000000000000e-07 + 0.00000000000000i,0.904262278204716 + 0.00363025663090923i;1.39000000000000e-07 + 0.00000000000000i,0.916566671040487 + 0.00352945275315903i;1.40000000000000e-07 + 0.00000000000000i,0.928581137328181 + 0.00343477596593924i;1.41000000000000e-07 + 0.00000000000000i,0.940315415996334 + 0.00334685618249554i;1.42000000000000e-07 + 0.00000000000000i,0.951786088128653 + 0.00326643226862334i;1.43000000000000e-07 + 0.00000000000000i,0.963001711542323 + 0.00319394661446027i;1.44000000000000e-07 + 0.00000000000000i,0.973964906981885 + 0.00312950045079895i;1.45000000000000e-07 + 0.00000000000000i,0.984698320022330 + 0.00307317783615491i;1.46000000000000e-07 + 0.00000000000000i,0.995201582399111 + 0.00302429820033668i;1.47000000000000e-07 + 0.00000000000000i,1.00548726758170 + 0.00298221933399098i;1.48000000000000e-07 + 0.00000000000000i,1.01556037764339 + 0.00294588051330091i;1.49000000000000e-07 + 0.00000000000000i,1.02543085563409 + 0.00291420227962127i;1.50000000000000e-07 + 0.00000000000000i,1.03510000333609 + 0.00288598948334354i;1.51000000000000e-07 + 0.00000000000000i,1.04457907486853 + 0.00286036616754394i;1.52000000000000e-07 + 0.00000000000000i,1.05387726583457 + 0.00283657563761663i;1.53000000000000e-07 + 0.00000000000000i,1.06299855846103 + 0.00281390724744562i;1.54000000000000e-07 + 0.00000000000000i,1.07194592588972 + 0.00279178403793741i;1.55000000000000e-07 + 0.00000000000000i,1.08073257270499 + 0.00277008183244213i;1.56000000000000e-07 + 0.00000000000000i,1.08935778067429 + 0.00274836999582396i;1.57000000000000e-07 + 0.00000000000000i,1.09783399131610 + 0.00272673208865188i;1.58000000000000e-07 + 0.00000000000000i,1.10616004450327 + 0.00270488718166149i;1.59000000000000e-07 + 0.00000000000000i,1.11434353560197 + 0.00268292610021559i;1.60000000000000e-07 + 0.00000000000000i,1.12238837461264 + 0.00266079428586223i;1.61000000000000e-07 + 0.00000000000000i,1.13030028397935 + 0.00263862031855042i;1.62000000000000e-07 + 0.00000000000000i,1.13807843209637 + 0.00261625475655988i;1.63000000000000e-07 + 0.00000000000000i,1.14573354873956 + 0.00259395858044178i;1.64000000000000e-07 + 0.00000000000000i,1.15325915335049 + 0.00257157397752373i;1.65000000000000e-07 + 0.00000000000000i,1.16067218880241 + 0.00254934400320423i;1.66000000000000e-07 + 0.00000000000000i,1.16797015739809 + 0.00252723109653152i;1.67000000000000e-07 + 0.00000000000000i,1.17515736668649 + 0.00250529570383340i;1.68000000000000e-07 + 0.00000000000000i,1.18223509061356 + 0.00248352737299786i;1.69000000000000e-07 + 0.00000000000000i,1.18920619430224 + 0.00246196571774666i;1.70000000000000e-07 + 0.00000000000000i,1.19607411551900 + 0.00244064612159652i;1.71000000000000e-07 + 0.00000000000000i,1.20284316653877 + 0.00241962646229877i;1.72000000000000e-07 + 0.00000000000000i,1.20951154778624 + 0.00239879907381242i;1.73000000000000e-07 + 0.00000000000000i,1.21608741729633 + 0.00237835756124418i;1.74000000000000e-07 + 0.00000000000000i,1.22256895645190 + 0.00235817136600961i;1.75000000000000e-07 + 0.00000000000000i,1.22896320624257 + 0.00233839867682677i;1.76000000000000e-07 + 0.00000000000000i,1.23526484586479 + 0.00231884162947911i;1.77000000000000e-07 + 0.00000000000000i,1.24148239869786 + 0.00229966391243556i;1.78000000000000e-07 + 0.00000000000000i,1.24761577957643 + 0.00228082534882069i;1.79000000000000e-07 + 0.00000000000000i,1.25366773335805 + 0.00226233547881962i;1.80000000000000e-07 + 0.00000000000000i,1.25964071023926 + 0.00224422057258406i;1.81000000000000e-07 + 0.00000000000000i,1.26553305501117 + 0.00222637946578820i;1.82000000000000e-07 + 0.00000000000000i,1.27135156044127 + 0.00220895374296009i;1.83000000000000e-07 + 0.00000000000000i,1.27709232997611 + 0.00219178704880326i;1.84000000000000e-07 + 0.00000000000000i,1.28276145916509 + 0.00217498925567409i;1.85000000000000e-07 + 0.00000000000000i,1.28835424641138 + 0.00215850582335836i;1.86000000000000e-07 + 0.00000000000000i,1.29388075421081 + 0.00214234317646779i;1.87000000000000e-07 + 0.00000000000000i,1.29934309815966 + 0.00212661420073837i;1.88000000000000e-07 + 0.00000000000000i,1.30473588437215 + 0.00211114527279883i;1.89000000000000e-07 + 0.00000000000000i,1.31006382720769 + 0.00209599799662808i;1.90000000000000e-07 + 0.00000000000000i,1.31532897049478 + 0.00208120330752432i;1.91000000000000e-07 + 0.00000000000000i,1.32052948508939 + 0.00206665209006173i;1.92000000000000e-07 + 0.00000000000000i,1.32566935148706 + 0.00205244374715600i;1.93000000000000e-07 + 0.00000000000000i,1.33075026500164 + 0.00203856197348573i;1.94000000000000e-07 + 0.00000000000000i,1.33577365027466 + 0.00202501462069634i;1.95000000000000e-07 + 0.00000000000000i,1.34073521373143 + 0.00201163123642667i;1.96000000000000e-07 + 0.00000000000000i,1.34564326822760 + 0.00199863016092670i;1.97000000000000e-07 + 0.00000000000000i,1.35049603931718 + 0.00198590648949382i;1.98000000000000e-07 + 0.00000000000000i,1.35529408730570 + 0.00197345156875452i;1.99000000000000e-07 + 0.00000000000000i,1.36004017544764 + 0.00196131417774930i;2.00000000000000e-07 + 0.00000000000000i,1.36473402541687 + 0.00194941449898592i;2.01000000000000e-07 + 0.00000000000000i,1.36937623938795 + 0.00193779201218026i;2.02000000000000e-07 + 0.00000000000000i,1.37396946574932 + 0.00192645013659301i;2.03000000000000e-07 + 0.00000000000000i,1.37851064791995 + 0.00191528982246805i;2.04000000000000e-07 + 0.00000000000000i,1.38299909789126 + 0.00190437803469285i;2.05000000000000e-07 + 0.00000000000000i,1.38744676167979 + 0.00189376041165146i;2.06000000000000e-07 + 0.00000000000000i,1.39184810113916 + 0.00188340753831382i;2.07000000000000e-07 + 0.00000000000000i,1.39620351850376 + 0.00187327973785263i;2.08000000000000e-07 + 0.00000000000000i,1.40051351013883 + 0.00186337215533602i;2.09000000000000e-07 + 0.00000000000000i,1.40477988141713 + 0.00185370228548127i;2.10000000000000e-07 + 0.00000000000000i,1.40900237088045 + 0.00184423768229084i;2.11000000000000e-07 + 0.00000000000000i,1.41318533107870 + 0.00183507467354655i;2.12000000000000e-07 + 0.00000000000000i,1.41732407557618 + 0.00182606047752684i;2.13000000000000e-07 + 0.00000000000000i,1.42142118548040 + 0.00181723672115282i;2.14000000000000e-07 + 0.00000000000000i,1.42547991259928 + 0.00180867891307771i;2.15000000000000e-07 + 0.00000000000000i,1.42949791722328 + 0.00180029775255298i;2.16000000000000e-07 + 0.00000000000000i,1.43347820168555 + 0.00179214935766731i;2.17000000000000e-07 + 0.00000000000000i,1.43741926458820 + 0.00178416096417493i;2.18000000000000e-07 + 0.00000000000000i,1.44132390973788 + 0.00177641890915369i;2.19000000000000e-07 + 0.00000000000000i,1.44519134929874 + 0.00176884574253726i;2.20000000000000e-07 + 0.00000000000000i,1.44902326929734 + 0.00176148465517056i;2.21000000000000e-07 + 0.00000000000000i,1.45281768975431 + 0.00175425768287853i;2.22000000000000e-07 + 0.00000000000000i,1.45657719256922 + 0.00174723489524373i;2.23000000000000e-07 + 0.00000000000000i,1.46030330253318 + 0.00174041581222902i;2.24000000000000e-07 + 0.00000000000000i,1.46399516719705 + 0.00173376146182461i;2.25000000000000e-07 + 0.00000000000000i,1.46764654252301 + 0.00172721846879393i;2.26000000000000e-07 + 0.00000000000000i,1.47127319967517 + 0.00172093579050301i;2.27000000000000e-07 + 0.00000000000000i,1.47486451935055 + 0.00171475505260568i;2.28000000000000e-07 + 0.00000000000000i,1.47842598108960 + 0.00170877852650954i;2.29000000000000e-07 + 0.00000000000000i,1.48195553831944 + 0.00170296256751984i;2.30000000000000e-07 + 0.00000000000000i,1.48545325711588 + 0.00169727516266127i;2.31000000000000e-07 + 0.00000000000000i,1.48892253726515 + 0.00169180253883112i;2.32000000000000e-07 + 0.00000000000000i,1.49236104958718 + 0.00168645567367402i;2.33000000000000e-07 + 0.00000000000000i,1.49576771954409 + 0.00168119811623351i;2.34000000000000e-07 + 0.00000000000000i,1.49914841240272 + 0.00167618665637289i;2.35000000000000e-07 + 0.00000000000000i,1.50250002904207 + 0.00167129118717089i;2.36000000000000e-07 + 0.00000000000000i,1.50582189369384 + 0.00166649702582776i;2.37000000000000e-07 + 0.00000000000000i,1.50911742005851 + 0.00166189117413798i;2.38000000000000e-07 + 0.00000000000000i,1.51238575230770 + 0.00165743273198301i;2.39000000000000e-07 + 0.00000000000000i,1.51562688260479 + 0.00165309971098296i;2.40000000000000e-07 + 0.00000000000000i,1.51884073706216 + 0.00164886825205636i;2.41000000000000e-07 + 0.00000000000000i,1.52202932322371 + 0.00164480294303481i;2.42000000000000e-07 + 0.00000000000000i,1.52519186700916 + 0.00164087707005275i;2.43000000000000e-07 + 0.00000000000000i,1.52832863692603 + 0.00163705767901107i;2.44000000000000e-07 + 0.00000000000000i,1.53143442910146 + 0.00163331410169481i;2.45000000000000e-07 + 0.00000000000000i,1.53452189318674 + 0.00162977018087066i;2.46000000000000e-07 + 0.00000000000000i,1.53758517128405 + 0.00162633489909321i;2.47000000000000e-07 + 0.00000000000000i,1.54062279017799 + 0.00162299500406882i;2.48000000000000e-07 + 0.00000000000000i,1.54363926168823 + 0.00161983165230185i;2.49000000000000e-07 + 0.00000000000000i,1.54663040274426 + 0.00161673707427661i;2.50000000000000e-07 + 0.00000000000000i,1.54959882101374 + 0.00161375743347656i;2.51000000000000e-07 + 0.00000000000000i,1.55254477578930 + 0.00161090842103939i;2.52000000000000e-07 + 0.00000000000000i,1.55546744055508 + 0.00160814599126118i;2.53000000000000e-07 + 0.00000000000000i,1.55836859055099 + 0.00160549816831682i;2.54000000000000e-07 + 0.00000000000000i,1.56124708083924 + 0.00160295364953226i;2.55000000000000e-07 + 0.00000000000000i,1.56410602608250 + 0.00160055968739388i;2.56000000000000e-07 + 0.00000000000000i,1.56694295617689 + 0.00159824902995772i;2.57000000000000e-07 + 0.00000000000000i,1.56975800019038 + 0.00159602329347550i;2.58000000000000e-07 + 0.00000000000000i,1.57255244382784 + 0.00159389280909877i;2.59000000000000e-07 + 0.00000000000000i,1.57532837972871 + 0.00159191393603527i;2.60000000000000e-07 + 0.00000000000000i,1.57808321478483 + 0.00159000504597359i;2.61000000000000e-07 + 0.00000000000000i,1.58081849550859 + 0.00158820836783098i;2.62000000000000e-07 + 0.00000000000000i,1.58353251508669 + 0.00158647544178586i;2.63000000000000e-07 + 0.00000000000000i,1.58622677087906 + 0.00158481302754717i;2.64000000000000e-07 + 0.00000000000000i,1.58890305502332 + 0.00158326970068982i;2.65000000000000e-07 + 0.00000000000000i,1.59155607435562 + 0.00158182373220053i;2.66000000000000e-07 + 0.00000000000000i,1.59419488327993 + 0.00158046502482762i;2.67000000000000e-07 + 0.00000000000000i,1.59681468567580 + 0.00157919061963943i;2.68000000000000e-07 + 0.00000000000000i,1.59941663581523 + 0.00157800380592777i;2.69000000000000e-07 + 0.00000000000000i,1.60200071742792 + 0.00157691038976884i;2.70000000000000e-07 + 0.00000000000000i,1.60456729020578 + 0.00157591060345532i;2.71000000000000e-07 + 0.00000000000000i,1.60711439958938 + 0.00157495024004416i;2.72000000000000e-07 + 0.00000000000000i,1.60964717069186 + 0.00157413570721475i;2.73000000000000e-07 + 0.00000000000000i,1.61215905289633 + 0.00157331812947024i;2.74000000000000e-07 + 0.00000000000000i,1.61465664058582 + 0.00157263926986708i;2.75000000000000e-07 + 0.00000000000000i,1.61713692533820 + 0.00157202431834975i;2.76000000000000e-07 + 0.00000000000000i,1.61960160716093 + 0.00157151288305887i;2.77000000000000e-07 + 0.00000000000000i,1.62204859531027 + 0.00157103984505773i;2.78000000000000e-07 + 0.00000000000000i,1.62448096109799 + 0.00157068837384148i;2.79000000000000e-07 + 0.00000000000000i,1.62689610610920 + 0.00157036610039431i;2.80000000000000e-07 + 0.00000000000000i,1.62929696336186 + 0.00157015992501630i;2.81000000000000e-07 + 0.00000000000000i,1.63168166305012 + 0.00157000146318895i;2.82000000000000e-07 + 0.00000000000000i,1.63404937731395 + 0.00156988841460586i;2.83000000000000e-07 + 0.00000000000000i,1.63640284765032 + 0.00156985961671164i;2.84000000000000e-07 + 0.00000000000000i,1.63873731986009 + 0.00156988982286948i;2.85000000000000e-07 + 0.00000000000000i,1.64106112296812 + 0.00157001247397319i;2.86000000000000e-07 + 0.00000000000000i,1.64336877860007 + 0.00157015828464525i;2.87000000000000e-07 + 0.00000000000000i,1.64566404226351 + 0.00157043741166694i;2.88000000000000e-07 + 0.00000000000000i,1.64794323686949 + 0.00157072448398582i;2.89000000000000e-07 + 0.00000000000000i,1.65020890092037 + 0.00157109737432713i;2.90000000000000e-07 + 0.00000000000000i,1.65246125148948 + 0.00157155138556034i;2.91000000000000e-07 + 0.00000000000000i,1.65469837154006 + 0.00157203043862454i;2.92000000000000e-07 + 0.00000000000000i,1.65692176139447 + 0.00157257130685859i;2.93000000000000e-07 + 0.00000000000000i,1.65913218387786 + 0.00157319401219555i;2.94000000000000e-07 + 0.00000000000000i,1.66132967355602 + 0.00157388603468852i;2.95000000000000e-07 + 0.00000000000000i,1.66351223940924 + 0.00157459380982688i;2.96000000000000e-07 + 0.00000000000000i,1.66568271876774 + 0.00157538676692760i;2.97000000000000e-07 + 0.00000000000000i,1.66784062854903 + 0.00157624887631064i;2.98000000000000e-07 + 0.00000000000000i,1.66998533654173 + 0.00157716395311941i;2.99000000000000e-07 + 0.00000000000000i,1.67211636981618 + 0.00157810140154477i;3.00000000000000e-07 + 0.00000000000000i,1.67423653141278 + 0.00157913562860027i;3.01000000000000e-07 + 0.00000000000000i,1.67634332092590 + 0.00158020640286474i;3.02000000000000e-07 + 0.00000000000000i,1.67843692028693 + 0.00158130944243515i;3.03000000000000e-07 + 0.00000000000000i,1.68051995825317 + 0.00158250995478999i;3.04000000000000e-07 + 0.00000000000000i,1.68258991014380 + 0.00158373469178891i;3.05000000000000e-07 + 0.00000000000000i,1.68464344119844 + 0.00158497440505099i;3.06000000000000e-07 + 0.00000000000000i,1.68668969803779 + 0.00158630693411371i;3.07000000000000e-07 + 0.00000000000000i,1.68872325481727 + 0.00158766611679433i;3.08000000000000e-07 + 0.00000000000000i,1.69074681547943 + 0.00158910782758728i;3.09000000000000e-07 + 0.00000000000000i,1.69275764058958 + 0.00159056795923741i;3.10000000000000e-07 + 0.00000000000000i,1.69475802462346 + 0.00159209280405929i;3.11000000000000e-07 + 0.00000000000000i,1.69674732848770 + 0.00159366915973188i;3.12000000000000e-07 + 0.00000000000000i,1.69872407575366 + 0.00159526195261920i;3.13000000000000e-07 + 0.00000000000000i,1.70069037857604 + 0.00159690746504201i;3.14000000000000e-07 + 0.00000000000000i,1.70264604200443 + 0.00159860924159604i;3.15000000000000e-07 + 0.00000000000000i,1.70459111467346 + 0.00160035951671720i;3.16000000000000e-07 + 0.00000000000000i,1.70652582113096 + 0.00160216148053175i;3.17000000000000e-07 + 0.00000000000000i,1.70844951203282 + 0.00160399678768548i;3.18000000000000e-07 + 0.00000000000000i,1.71036203397876 + 0.00160585306936378i;3.19000000000000e-07 + 0.00000000000000i,1.71226463493629 + 0.00160776310206718i;3.20000000000000e-07 + 0.00000000000000i,1.71415931149796 + 0.00160977611942142i;3.21000000000000e-07 + 0.00000000000000i,1.71604182157145 + 0.00161177901872948i;3.22000000000000e-07 + 0.00000000000000i,1.71791386841116 + 0.00161381363430948i;3.23000000000000e-07 + 0.00000000000000i,1.71977604255343 + 0.00161588821602338i;3.24000000000000e-07 + 0.00000000000000i,1.72162860093711 + 0.00161801203277580i;3.25000000000000e-07 + 0.00000000000000i,1.72346708435759 + 0.00162015382661928i;3.26000000000000e-07 + 0.00000000000000i,1.72529928498375 + 0.00162233540130808i;3.27000000000000e-07 + 0.00000000000000i,1.72712334764071 + 0.00162458395517568i;3.28000000000000e-07 + 0.00000000000000i,1.72893742908816 + 0.00162686152159581i;3.29000000000000e-07 + 0.00000000000000i,1.73074191180434 + 0.00162917628488905i;3.30000000000000e-07 + 0.00000000000000i,1.73253747471961 + 0.00163153742289443i;3.31000000000000e-07 + 0.00000000000000i,1.73432360900534 + 0.00163392897134024i;3.32000000000000e-07 + 0.00000000000000i,1.73609985939283 + 0.00163635113558378i;3.33000000000000e-07 + 0.00000000000000i,1.73786715088572 + 0.00163879271782394i;3.34000000000000e-07 + 0.00000000000000i,1.73962653623745 + 0.00164129176085048i;3.35000000000000e-07 + 0.00000000000000i,1.74137671751901 + 0.00164382890025260i;3.36000000000000e-07 + 0.00000000000000i,1.74311817670286 + 0.00164640105427221i;3.37000000000000e-07 + 0.00000000000000i,1.74485121946114 + 0.00164901392638813i;3.38000000000000e-07 + 0.00000000000000i,1.74657458289201 + 0.00165163064563602i;3.39000000000000e-07 + 0.00000000000000i,1.74828992399855 + 0.00165429649951866i;3.40000000000000e-07 + 0.00000000000000i,1.74999795329406 + 0.00165703019182971i;3.41000000000000e-07 + 0.00000000000000i,1.75169618285177 + 0.00165974943215196i;3.42000000000000e-07 + 0.00000000000000i,1.75338654006288 + 0.00166252170172734i;3.43000000000000e-07 + 0.00000000000000i,1.75506892950072 + 0.00166533404110521i;3.44000000000000e-07 + 0.00000000000000i,1.75674315220322 + 0.00166817570235844i;3.45000000000000e-07 + 0.00000000000000i,1.75840504569909 + 0.00167101652245911i;3.46000000000000e-07 + 0.00000000000000i,1.76006265594369 + 0.00167390508165124i;3.47000000000000e-07 + 0.00000000000000i,1.76171217625986 + 0.00167682690320242i;3.48000000000000e-07 + 0.00000000000000i,1.76335487859725 + 0.00167980650764906i;3.49000000000000e-07 + 0.00000000000000i,1.76498789691819 + 0.00168276644181965i;3.50000000000000e-07 + 0.00000000000000i,1.76661499931462 + 0.00168580331996920i;3.51000000000000e-07 + 0.00000000000000i,1.76823460605867 + 0.00168887072906554i;3.52000000000000e-07 + 0.00000000000000i,1.76984607566272 + 0.00169195467239880i;3.53000000000000e-07 + 0.00000000000000i,1.77144949356476 + 0.00169506070234602i;3.54000000000000e-07 + 0.00000000000000i,1.77304600307999 + 0.00169819542342595i;3.55000000000000e-07 + 0.00000000000000i,1.77463474882252 + 0.00170136113836678i;3.56000000000000e-07 + 0.00000000000000i,1.77621706636296 + 0.00170458039230982i;3.57000000000000e-07 + 0.00000000000000i,1.77779129030194 + 0.00170779964054975i;3.58000000000000e-07 + 0.00000000000000i,1.77935954945264 + 0.00171107580737343i;3.59000000000000e-07 + 0.00000000000000i,1.78091964040242 + 0.00171435460896266i;3.60000000000000e-07 + 0.00000000000000i,1.78247423697468 + 0.00171770527073197i;3.61000000000000e-07 + 0.00000000000000i,1.78402052214764 + 0.00172104671776108i;3.62000000000000e-07 + 0.00000000000000i,1.78555954131689 + 0.00172441360721107i;3.63000000000000e-07 + 0.00000000000000i,1.78709293569200 + 0.00172783800870560i;3.64000000000000e-07 + 0.00000000000000i,1.78861894767262 + 0.00173127201461156i;3.65000000000000e-07 + 0.00000000000000i,1.79013457615075 + 0.00173472020015054i;3.66000000000000e-07 + 0.00000000000000i,1.79164602857671 + 0.00173817745260830i;3.67000000000000e-07 + 0.00000000000000i,1.79315287316677 + 0.00174171921078980i;3.68000000000000e-07 + 0.00000000000000i,1.79465158539083 + 0.00174524289740877i;3.69000000000000e-07 + 0.00000000000000i,1.79614533928393 + 0.00174883853353015i;3.70000000000000e-07 + 0.00000000000000i,1.79763098109059 + 0.00175241551753188i;3.71000000000000e-07 + 0.00000000000000i,1.79911136135302 + 0.00175604873360497i;3.72000000000000e-07 + 0.00000000000000i,1.80058542092365 + 0.00175970861180865i;3.73000000000000e-07 + 0.00000000000000i,1.80205356479602 + 0.00176340345138778i;3.74000000000000e-07 + 0.00000000000000i,1.80351483463778 + 0.00176710879987267i;3.75000000000000e-07 + 0.00000000000000i,1.80497009811034 + 0.00177084330661860i;3.76000000000000e-07 + 0.00000000000000i,1.80641927782992 + 0.00177461055571334i;3.77000000000000e-07 + 0.00000000000000i,1.80786335723693 + 0.00177841677897691i;3.78000000000000e-07 + 0.00000000000000i,1.80930069450027 + 0.00178223693706174i;3.79000000000000e-07 + 0.00000000000000i,1.81073096050082 + 0.00178605646608192i;3.80000000000000e-07 + 0.00000000000000i,1.81215865822561 + 0.00178997809532586i;3.81000000000000e-07 + 0.00000000000000i,1.81357780135212 + 0.00179385761074358i;3.82000000000000e-07 + 0.00000000000000i,1.81499141565132 + 0.00179777892677088i;3.83000000000000e-07 + 0.00000000000000i,1.81640052334325 + 0.00180175999300615i;3.84000000000000e-07 + 0.00000000000000i,1.81780402404335 + 0.00180575021168704i;3.85000000000000e-07 + 0.00000000000000i,1.81919760120533 + 0.00180974501768114i;3.86000000000000e-07 + 0.00000000000000i,1.82058822719985 + 0.00181376542198423i;3.87000000000000e-07 + 0.00000000000000i,1.82197434388103 + 0.00181782591457580i;3.88000000000000e-07 + 0.00000000000000i,1.82335557136247 + 0.00182193209207795i;3.89000000000000e-07 + 0.00000000000000i,1.82473069646856 + 0.00182604920843830i;3.90000000000000e-07 + 0.00000000000000i,1.82610037551661 + 0.00183018857846406i;3.91000000000000e-07 + 0.00000000000000i,1.82746499105908 + 0.00183436642470471i;3.92000000000000e-07 + 0.00000000000000i,1.82882440436611 + 0.00183856485499995i;3.93000000000000e-07 + 0.00000000000000i,1.83017900239235 + 0.00184280070234464i;3.94000000000000e-07 + 0.00000000000000i,1.83152783951891 + 0.00184704843195760i;3.95000000000000e-07 + 0.00000000000000i,1.83287166091667 + 0.00185133629449803i;3.96000000000000e-07 + 0.00000000000000i,1.83421175930944 + 0.00185567079186318i;3.97000000000000e-07 + 0.00000000000000i,1.83554477431921 + 0.00185999192747463i;3.98000000000000e-07 + 0.00000000000000i,1.83687384628808 + 0.00186435407361837i;3.99000000000000e-07 + 0.00000000000000i,1.83819835053717 + 0.00186876382083647i;4.00000000000000e-07 + 0.00000000000000i,1.83951813938881 + 0.00187318950480814i;4.01000000000000e-07 + 0.00000000000000i,1.84083223923085 + 0.00187762654266198i;4.02000000000000e-07 + 0.00000000000000i,1.84214192269206 + 0.00188211383214196i;4.03000000000000e-07 + 0.00000000000000i,1.84344667851421 + 0.00188662249442907i;4.04000000000000e-07 + 0.00000000000000i,1.84474768346254 + 0.00189116764666866i;4.05000000000000e-07 + 0.00000000000000i,1.84603934094339 + 0.00189571044274829i;4.06000000000000e-07 + 0.00000000000000i,1.84732961192391 + 0.00190028614600904i;4.07000000000000e-07 + 0.00000000000000i,1.84861635589642 + 0.00190491194710596i;4.08000000000000e-07 + 0.00000000000000i,1.84989864003960 + 0.00190956806861091i;4.09000000000000e-07 + 0.00000000000000i,1.85117548650525 + 0.00191424282976723i;4.10000000000000e-07 + 0.00000000000000i,1.85244804421984 + 0.00191893903534342i;4.11000000000000e-07 + 0.00000000000000i,1.85371644161474 + 0.00192367156498006i;4.12000000000000e-07 + 0.00000000000000i,1.85498069234043 + 0.00192843928082091i;4.13000000000000e-07 + 0.00000000000000i,1.85624051946303 + 0.00193322662446963i;4.14000000000000e-07 + 0.00000000000000i,1.85749593505771 + 0.00193804772621710i;4.15000000000000e-07 + 0.00000000000000i,1.85874801042411 + 0.00194291918440149i;4.16000000000000e-07 + 0.00000000000000i,1.85999478392622 + 0.00194779012807034i;4.17000000000000e-07 + 0.00000000000000i,1.86123743765114 + 0.00195269232613460i;4.18000000000000e-07 + 0.00000000000000i,1.86247631628756 + 0.00195762799138492i;4.19000000000000e-07 + 0.00000000000000i,1.86371143989375 + 0.00196260536476470i;4.20000000000000e-07 + 0.00000000000000i,1.86494320154036 + 0.00196761943006954i;4.21000000000000e-07 + 0.00000000000000i,1.86617051068579 + 0.00197266963925333i;4.22000000000000e-07 + 0.00000000000000i,1.86739245083495 + 0.00197770775791918i;4.23000000000000e-07 + 0.00000000000000i,1.86861139353649 + 0.00198280621838418i;4.24000000000000e-07 + 0.00000000000000i,1.86982580804323 + 0.00198791018856766i;4.25000000000000e-07 + 0.00000000000000i,1.87103466849501 + 0.00199306771068553i;4.26000000000000e-07 + 0.00000000000000i,1.87224125542170 + 0.00199823887320934i;4.27000000000000e-07 + 0.00000000000000i,1.87344565796057 + 0.00200346878398305i;4.28000000000000e-07 + 0.00000000000000i,1.87464445492667 + 0.00200869167232515i;4.29000000000000e-07 + 0.00000000000000i,1.87584059541312 + 0.00201398814495861i;4.30000000000000e-07 + 0.00000000000000i,1.87703322879025 + 0.00201929475711467i;4.31000000000000e-07 + 0.00000000000000i,1.87822196325483 + 0.00202463854975876i;4.32000000000000e-07 + 0.00000000000000i,1.87940638663708 + 0.00202999445117827i;4.33000000000000e-07 + 0.00000000000000i,1.88058910075927 + 0.00203541969307725i;4.34000000000000e-07 + 0.00000000000000i,1.88176704478305 + 0.00204085893606199i;4.35000000000000e-07 + 0.00000000000000i,1.88294129365089 + 0.00204632601365404i;4.36000000000000e-07 + 0.00000000000000i,1.88411289425242 + 0.00205184263256265i;4.37000000000000e-07 + 0.00000000000000i,1.88528036057668 + 0.00205736173659018i;4.38000000000000e-07 + 0.00000000000000i,1.88644488342914 + 0.00206293932771733i;4.39000000000000e-07 + 0.00000000000000i,1.88760532040701 + 0.00206853029796271i;4.40000000000000e-07 + 0.00000000000000i,1.88876428960351 + 0.00207418767478174i;4.41000000000000e-07 + 0.00000000000000i,1.88991825227024 + 0.00207984790857667i;4.42000000000000e-07 + 0.00000000000000i,1.89106930845336 + 0.00208555721727697i;4.43000000000000e-07 + 0.00000000000000i,1.89221701658997 + 0.00209129074650746i;4.44000000000000e-07 + 0.00000000000000i,1.89336215650669 + 0.00209707389642466i;4.45000000000000e-07 + 0.00000000000000i,1.89449910578987 + 0.00210283139142779i;4.46000000000000e-07 + 0.00000000000000i,1.89563818968538 + 0.00210868925586485i;4.47000000000000e-07 + 0.00000000000000i,1.89677301710307 + 0.00211454709254603i;4.48000000000000e-07 + 0.00000000000000i,1.89790483768165 + 0.00212045736567504i;4.49000000000000e-07 + 0.00000000000000i,1.89903343816309 + 0.00212639502458929i;4.50000000000000e-07 + 0.00000000000000i,1.90015915837079 + 0.00213236287580546i;4.51000000000000e-07 + 0.00000000000000i,1.90128247330844 + 0.00213838163079079i;4.52000000000000e-07 + 0.00000000000000i,1.90240309345652 + 0.00214444255439813i;4.53000000000000e-07 + 0.00000000000000i,1.90352015415906 + 0.00215052654486842i;4.54000000000000e-07 + 0.00000000000000i,1.90463444452861 + 0.00215664695042360i;4.55000000000000e-07 + 0.00000000000000i,1.90574575036597 + 0.00216280249571497i;4.56000000000000e-07 + 0.00000000000000i,1.90685410757933 + 0.00216898461677267i;4.57000000000000e-07 + 0.00000000000000i,1.90795970828038 + 0.00217521243856006i;4.58000000000000e-07 + 0.00000000000000i,1.90906210072259 + 0.00218146876774768i;4.59000000000000e-07 + 0.00000000000000i,1.91016279221799 + 0.00218777395682846i;4.60000000000000e-07 + 0.00000000000000i,1.91126142827220 + 0.00219413536513543i;4.61000000000000e-07 + 0.00000000000000i,1.91235604490826 + 0.00220050513937833i;4.62000000000000e-07 + 0.00000000000000i,1.91344752535900 + 0.00220690214953608i;4.63000000000000e-07 + 0.00000000000000i,1.91453619332871 + 0.00221333889310210i;4.64000000000000e-07 + 0.00000000000000i,1.91562445326088 + 0.00221984742675903i;4.65000000000000e-07 + 0.00000000000000i,1.91670467400591 + 0.00222634556562395i;4.66000000000000e-07 + 0.00000000000000i,1.91778610455611 + 0.00223289553059637i;4.67000000000000e-07 + 0.00000000000000i,1.91886506776688 + 0.00223949458739107i;4.68000000000000e-07 + 0.00000000000000i,1.91994150004208 + 0.00224613467628764i;4.69000000000000e-07 + 0.00000000000000i,1.92101570091332 + 0.00225280716149084i;4.70000000000000e-07 + 0.00000000000000i,1.92208692846140 + 0.00225951747809027i;4.71000000000000e-07 + 0.00000000000000i,1.92315635001696 + 0.00226626520288901i;4.72000000000000e-07 + 0.00000000000000i,1.92422351402071 + 0.00227306532987394i;4.73000000000000e-07 + 0.00000000000000i,1.92528762716518 + 0.00227989894040512i;4.74000000000000e-07 + 0.00000000000000i,1.92634934433526 + 0.00228675646141861i;4.75000000000000e-07 + 0.00000000000000i,1.92740859900941 + 0.00229364939496724i;4.76000000000000e-07 + 0.00000000000000i,1.92846677895890 + 0.00230061763771818i;4.77000000000000e-07 + 0.00000000000000i,1.92952125015339 + 0.00230758378947711i;4.78000000000000e-07 + 0.00000000000000i,1.93057405163132 + 0.00231460685604070i;4.79000000000000e-07 + 0.00000000000000i,1.93162542297991 + 0.00232169105945000i;4.80000000000000e-07 + 0.00000000000000i,1.93267370711246 + 0.00232879353746997i;4.81000000000000e-07 + 0.00000000000000i,1.93372059790029 + 0.00233594590684591i;4.82000000000000e-07 + 0.00000000000000i,1.93476405201482 + 0.00234312005103159i;4.83000000000000e-07 + 0.00000000000000i,1.93580485899771 + 0.00235032459950016i;4.84000000000000e-07 + 0.00000000000000i,1.93684423546494 + 0.00235758675797081i;4.85000000000000e-07 + 0.00000000000000i,1.93787883433320 + 0.00236487484967514i;4.86000000000000e-07 + 0.00000000000000i,1.93891489653181 + 0.00237222737655468i;4.87000000000000e-07 + 0.00000000000000i,1.93994738852163 + 0.00237958950158838i;4.88000000000000e-07 + 0.00000000000000i,1.94097862549330 + 0.00238701379964723i;4.89000000000000e-07 + 0.00000000000000i,1.94200718887187 + 0.00239446561611091i;4.90000000000000e-07 + 0.00000000000000i,1.94303472398668 + 0.00240197961772917i;4.91000000000000e-07 + 0.00000000000000i,1.94406038985340 + 0.00240952841446600i;4.92000000000000e-07 + 0.00000000000000i,1.94508302509459 + 0.00241709267519817i;4.93000000000000e-07 + 0.00000000000000i,1.94610398945785 + 0.00242471840161723i;4.94000000000000e-07 + 0.00000000000000i,1.94712384434942 + 0.00243239337621903i;4.95000000000000e-07 + 0.00000000000000i,1.94814117118097 + 0.00244009755761966i;4.96000000000000e-07 + 0.00000000000000i,1.94915713636340 + 0.00244785067999597i;4.97000000000000e-07 + 0.00000000000000i,1.95017058790609 + 0.00245563176400173i;4.98000000000000e-07 + 0.00000000000000i,1.95118249396880 + 0.00246346305968352i;4.99000000000000e-07 + 0.00000000000000i,1.95219308317003 + 0.00247133821167010i;5.00000000000000e-07 + 0.00000000000000i,1.95320155221684 + 0.00247925516383704i];
d = Data(:,1);neff = (Data(:,2));perm = neff.^2;

figure;plot(d*1e9,real(perm));
figure;plot(d*1e9,imag(perm))

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
legend(["Imaginary permittivity","Real permittivity"])
title('PECVD 90nm SiN');

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

figure;yyaxis left;
plot(z*1e6,real(eps));set(gcf, 'Position', [00, 00, 400, 300]);xlabel('x(\mum)');ylabel('Re[\epsilon]');yyaxis right
plot(z*1e6,imag(eps));set(gcf, 'Position', [00, 00, 400, 300]);xlabel('x(\mum)');ylabel('Im[\epsilon]');

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
SiO_Y = [d_Al + 180e-9 -d_Al - 180e-9]/2*1e6;
SiO_Z = [z -z];


% % %%%%3.3^2 +(0.03*exp(-((z-5.168[um])/2[um])^2) + 0.03*exp(-((z+5.168[um])/2[um])^2) + 0.017*exp(-(z/10[um])^2))*i
% % tm = 0.03*exp(-((z-5.168e-6)/2e-6).^2) + 0.03*exp(-((z+5.168e-6)/2e-6).^2) + 0.017*exp(-(z/10e-6).^2);
% % figure
% % plot(z,tm)