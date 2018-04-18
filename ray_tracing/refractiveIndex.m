function n = refractiveIndex(material, x)
%All formulas from refractiveindex.info unless otherwise specified.
%Wavelength, X, should be in nanometers

if iscell(material)
    n(length(material),length(x)) = 0;
    m = material;
else
    n(length(x)) = 0;
    m{1} = material;
end

x = x*0.001;
l = size(n,1);


for i = 1:l

    switch m{i}
        case 'N-BK7'
            n(i,:) = sqrt( 1 + 1.03961212*x.^2./(x.^2-0.00600069867) + 0.231792344*x.^2./(x.^2-0.0200179144) + 1.01046945*x.^2./(x.^2-103.560653) );
        case 'FusedSilica'
            n(i,:) = sqrt( 1 + 0.6961663*x.^2./(x.^2-0.0684043^2) + 0.4079426*x.^2./(x.^2-0.1162414^2) + 0.8974794*x.^2./(x.^2-9.896161^2) );
        case 'SF2'
            n(i,:) = sqrt( 1 + 1.40301821*x.^2./(x.^2-0.0105795466) + 0.231767504*x.^2./(x.^2-0.0493226978) + 0.939056586*x.^2./(x.^2-112.405955) );
        case 'SF5'
            n(i,:) = sqrt( 1 + 1.46141885*x.^2./(x.^2-0.0111826126) + 0.247713019*x.^2./(x.^2-0.0508594669) + 0.949995832*x.^2./(x.^2-112.041888) );
        case 'SF10'
            n(i,:) = sqrt( 1 + 1.61625977*x.^2./(x.^2-0.0127534559) + 0.259229334*x.^2./(x.^2-0.0581983954) + 1.07762317*x.^2./(x.^2-116.60768) );
        case 'N-SF5'
            n(i,:) = sqrt( 1 + 1.52481889*x.^2./(x.^2-0.011254756) + 0.187085527*x.^2./(x.^2-0.0588995392) + 1.42729015*x.^2./(x.^2-129.141675) );
        case 'N-SF6HT'
            n(i,:) = sqrt( 1 + 1.77931763*x.^2./(x.^2-0.0133714182) + 0.338149866*x.^2./(x.^2-0.0617533621) + 2.08734474*x.^2./(x.^2-174.01759) );
        case 'N-SF11'
            n(i,:) = sqrt( 1 + 1.73759695*x.^2./(x.^2-0.013188707) + 0.313747346*x.^2./(x.^2-0.0623068142) + 1.89878101*x.^2./(x.^2-155.23629) );
        case 'N-BAF10'
            n(i,:) = sqrt( 1 + 1.5851495*x.^2./(x.^2-0.00926681282) + 0.143559385*x.^2./(x.^2-0.0424489805) + 1.08521269*x.^2./(x.^2-105.613573) );
        case 'E-BAF11'
            n(i,:) = sqrt( 2.71954649 + -1.00472501E-2*x.^2 + 2.00301385E-2*x.^(-2) + 4.65868302E-4*x.^(-4) + -7.51633336E-6*x.^(-6) + 1.77544989E-6*x.^(-8) );
        case 'N-BAK4'
            n(i,:) = sqrt( 1 + 1.28834642*x.^2./(x.^2-0.00779980626) + 0.132817724*x.^2./(x.^2-0.0315631177) + 0.945395373*x.^2./(x.^2-105.965875) );
        case 'FD10' %E-FD10
            n(i,:) = sqrt( 2.8815180 + -1.3228312E-2*x.^2 + 3.1455590E-2*x.^(-2) + 2.6851666E-3*x.^(-4) + -2.2577544E-4*x.^(-6) + 2.4693268E-5*x.^(-8) );
        case 'N-SSK5'
            n(i,:) = sqrt( 1 + 1.59222659*x.^2./(x.^2-0.00920284626) + 0.103520774*x.^2./(x.^2-0.0423530072) + 1.05174016*x.^2./(x.^2-106.927374) );
        case 'LAFN7'
            n(i,:) = sqrt( 1 + 1.66842615*x.^2./(x.^2-0.0103159999) + 0.298512803*x.^2./(x.^2-0.0469216348) + 1.0774376*x.^2./(x.^2-82.5078509) );
        case 'B270'
            % http://www.crystran.co.uk/userfiles/files/optical-glass-b270-refractive-index.pdf
            Availible_n = [706.5 1.51883;
                656.3 1.52037;
                643.8 1.52080;
                632.8 1.51120;
                589.3 1.52299;
                587.6 1.52307;
                546.1 1.52520;
                486.1 1.52929;
                480.0 1.52980;
                435.8 1.53416;
                404.7 1.53820;
                365.0 1.54510];
            x = x*1000;
            
            minIdx = x < Availible_n(end,1);
            maxIdx = x > Availible_n(1,1);
            n(i,minIdx) = Availible_n(end,2);
            n(i,maxIdx) = Availible_n(1,2);
            n(i,~minIdx & ~maxIdx) = interp1(Availible_n(:,1),Availible_n(:,2),x(~minIdx & ~maxIdx));
%             idxMin = x<Availible_n(:,1);
%             exact = find(Availible_n(:,1)-x==0,1);
%             if sum(idxMin)==0 %if x is greater than all of the wavelengths assign the index for the greatest wavelength
%                 n(i,:) = Availible_n(1,2);
%             elseif all(idxMin) %if x is smaller than all of the wavelengths, assign the index for the smallest wavelength
%                 n(i,:) = Availible_n(end,2);
%             elseif ~isempty(exact)
%                 n(i,:) = Availible_n(exact,2);
%             else
%                 idx = find(idxMin,1,'last');
%                 ratio = (Availible_n(idx,1)-x)/(Availible_n(idx,1)-Availible_n(idx+1,1));
%                 n(i,:) = Availible_n(idx,2) - ratio*(Availible_n(idx,2)-Availible_n(idx+1,2));
%             end
        case {'Aluminum', 'Al'}
            % Data set Rakic 1998: n,k 0.2066-12.4 um
            dat = [0.2066 0.122846740621+2.2878524992i;0.21089508319 0.12766323937+2.34380200581i;0.215279458439 0.132638075903+2.40069491191i;0.219754982075 0.137780447014+2.45856466247i;0.224323549012 0.143101094237+2.51744498718i;0.228987093564 0.148612342484+2.57736980208i;0.233747590255 0.154328130258+2.63837311416i;0.238607054657 0.160264034101+2.70048892948i;0.243567544248 0.166437289906+2.76375116492i;0.248631159276 0.172866813664+2.82819356301i;0.253800043655 0.179573223987+2.8938496086i;0.259076385867 0.186578868465+2.96075244558i;0.264462419893 0.193907855581+3.0289347913i;0.269960426156 0.201586093411+3.09842884581i;0.275572732488 0.209641335877+3.16926619271i;0.281301715116 0.218103236684+3.2414776876i;0.287149799665 0.227003410362+3.31509333009i;0.293119462189 0.23637549904+3.3901421147i;0.299213230217 0.246255242528+3.46665185552i;0.305433683824 0.256680548133+3.54464897941i;0.311783456723 0.267691555141+3.62415828216i;0.318265237381 0.279330687127+3.70520264187i;0.324881770154 0.291642683088+3.78780268424i;0.331635856454 0.304674595797+3.87197639482i;0.338530355932 0.318475742693+3.9577386747i;0.345568187689 0.333097591082+4.04510083794i;0.352752331512 0.348593555506+4.13407005265i;0.360085829136 0.36501868105+4.2246487324i;0.367571785533 0.38242918251+4.31683389164i;0.375213370223 0.40088180634+4.41061648852i;0.383013818621 0.420432981139+4.50598079085i;0.390976433402 0.441137724284+4.60290381618i;0.399104585903 0.463048278801+4.70135491513i;0.407401717548 0.486212467254+4.80129558656i;0.415871341306 0.510671770022+4.90267963255i;0.42451704318 0.536459164782+5.00545377803i;0.433342483722 0.563596802793+5.10955889075i;0.442351399584 0.592093644475+5.21493193877i;0.451547605104 0.621943229657+5.32150881124i;0.460934993913 0.653121812831+5.4292281021i;0.470517540592 0.685587146407+5.53803591403i;0.480299302349 0.719278240973+5.64789168451i;0.49028442074 0.754116467921+5.75877496865i;0.50047712342 0.790008396569+5.87069304304i;0.510881725935 0.826850779326+5.98368912413i;0.521502633549 0.864538124501+6.09785093007i;0.532344343107 0.902973343951+6.21331925589i;0.543411444943 0.942082056657+6.33029617015i;0.554708624818 0.981831304254+6.4490523562i;0.56624066591 1.02225373683+6.56993296731i;0.578012450835 1.06347881745+6.69336105689i;0.590028963715 1.1057733467+6.81983702415i;0.602295292289 1.14959470086+6.94993128614i;0.614816630068 1.19566164222+7.08426499818i;0.627598278532 1.24504922419+7.22346907649i;0.640645649374 1.29931534855+7.36810319041i;0.653964266795 1.36066421788+7.51850073134i;0.667559769839 1.43213865299+7.67447901036i;0.681437914784 1.51779027773+7.83481498635i;0.695604577575 1.62266619156+7.99635410453i;0.710065756315 1.75222713279+8.15267647823i;0.724827573806 1.91052144009+8.29259788304i;0.739896280135 2.09647891709+8.39965903025i;0.755278255329 2.29890388833+8.4547767539i;0.770980012047 2.49334714318+8.44368443416i;0.787008198346 2.64561743502+8.36704341358i;0.803369600489 2.72349127944+8.24645006553i;0.820071145824 2.71103705761+8.11992069579i;0.837119905712 2.61596564261+8.02772059943i;0.854523098523 2.46478082946+7.99728212025i;0.872288092693 2.28972164363+8.03655549328i;0.890422409845 2.11679793746+8.13783713444i;0.908933727968 1.96100825745+8.28647676114i;0.927829884674 1.8278962783+8.46796644208i;0.947118880514 1.71735276624+8.67111122492i;0.966808882365 1.62680472747+8.88837171746i;0.986908226887 1.55306618966+9.11507096386i;1.00742542406 1.49316794145+9.34846325893i;1.02836916076 1.44462731533+9.58700960496i;1.0497483045 1.40546877376+9.82989388106i;1.07157190709 1.37415565333+10.0767222422i;1.09384920857 1.34950329164+10.3273430211i;1.11658964103 1.33059968936+10.5817408634i;1.13980283269 1.31674085707+10.8399752049i;1.16349861191 1.3073807396+11.1021448578i;1.18768701137 1.30209335641+11.3683678742i;1.21237827235 1.30054448596+11.6387703118i;1.237582849 1.30247056289+11.9134801698i;1.26331141285 1.30766293993+12.1926243006i;1.28957485725 1.31595611267+12.4763270151i;1.31638430203 1.32721886693+12.7647096249i;1.34375109819 1.34134758426+13.0578904823i;1.37168683271 1.35826114602+13.3559852631i;1.40020333347 1.37789702659+13.6591073453i;1.42931267422 1.40020827547+13.9673682052i;1.45902717974 1.42516116713+14.2808777875i;1.48935943101 1.45273335536+14.59974483i;1.52032227059 1.48291241079+14.9240771363i;1.55192880799 1.51569465093+15.2539817956i;1.58419242529 1.55108419495+15.5895653531i;1.61712678277 1.58909219173+15.9309339374i;1.65074582469 1.62973618228+16.2781933486i;1.6850637852 1.6730395668+16.6314491137i;1.72009519439 1.71903115353+16.9908065152i;1.75585488439 1.76774477148+17.356370595i;1.7923579957 1.81921893348+17.7282461402i;1.82961998359 1.87349653834+18.1065376524i;1.86765662461 1.93062460372+18.4913493038i;1.90648402331 1.99065402252+18.8827848838i;1.94611861905 2.0536393373+19.2809477359i;1.98657719294 2.11963852802+19.6859406885i;2.02787687497 2.18871280917+20.0978659816i;2.07003515123 2.26092643313+20.5168251881i;2.11306987137 2.33634649686+20.9429191351i;2.15699925609 2.41504274957+21.3762478222i;2.2018419049 2.49708739927+21.8169103416i;2.24761680399 2.58255491623+22.2650047988i;2.29434333425 2.67152183188+22.7206282364i;2.34204127949 2.76406653149+23.1838765617i;2.39073083481 2.86026903944+23.6548444798i;2.44043261516 2.96021079591+24.1336254333i;2.49116766405 3.06397442393+24.6203115503i;2.54295746248 3.17164348604+25.114993603i;2.59582393804 3.28330222977+25.6177609772i;2.64978947414 3.39903532142+26.1287016561i;2.70487691958 3.51892756799+26.6479022196i;2.76110959812 3.64306362677+27.1754478621i;2.81851131845 3.77152770315+27.7114224299i;2.8771063842 3.90440323655+28.2559084822i;2.93691960427 4.04177257547+28.8089873774i;2.99797630331 4.18371664239+29.3707393872i;3.06030233246 4.33031458995+29.9412438429i;3.1239240803 4.48164345005+30.5205793145i;3.188868484 4.63777777801+31.1088238278i;3.25516304072 4.79878929429+31.7060551214i;3.32283581931 4.96474652702+32.3123509463i;3.39191547211 5.13571445874+32.9277894117i;3.46243124716 5.3117541816+33.5524493798i;3.53441300053 5.49292256588+34.1864109106i;3.60789120897 5.67927194695+34.829755761i;3.68289698284 5.87084983687+35.4825679384i;3.75946207926 6.06769866707+36.1449343106i;3.83761891556 6.26985556939+36.8169452721i;3.91740058299 6.4773522031+37.4986954682i;3.99884086078 6.6902146363+38.1902845733i;4.08197423038 6.90846329018+38.8918181236i;4.1668358901 7.13211295542+39.6034084007i;4.25346176999 7.36117288966+40.3251753601i;4.34188854708 7.59564700545+41.0572476008i;4.43215366088 7.83553415787+41.7997633679i;4.52429532923 8.08082854042+42.5528715791i;4.61835256454 8.33152019749+43.3167328656i;4.71436519021 8.58759566082+44.0915206146i;4.81237385758 8.84903871599+44.8774220006i;4.91242006309 9.11583130396+45.6746389906i;5.01454616587 9.38795456074+46.4833893057i;5.11879540566 9.66538999618+47.3039073236i;5.22521192114 9.94812081105+48.1364449016i;5.3338407686 10.2361333486+48.9812721034i;5.44472794102 10.5294186743+49.8386778074i;5.55792038754 10.8279742751+50.7089701815i;5.67346603336 11.1318058656+51.5924770013i;5.79141379999 11.4409292866+52.489545799i;5.91181362602 11.755372479+53.400543823i;6.0347164882 12.0751775109+54.3258577976i;6.16017442306 12.4004026381+55.2658934693i;6.28824054896 12.7311243705+56.2210749328i;6.41896908852 13.0674395201+57.1918437304i;6.55241539166 13.4094672024+58.1786577242i;6.68863595894 13.7573507641+59.1819897428i;6.82768846557 14.1112596082+60.2023260089i;6.96963178576 14.4713908912+61.2401643569i;7.11452601772 14.837971065+62.296012256i;7.26243250901 15.2112572405+63.3703846539i;7.41341388261 15.5915383519+64.4638016636i;7.56753406337 15.9791361024+65.576786115i;7.72485830511 16.3744056768+66.7098609976i;7.88545321822 16.7777362097+67.8635468221i;8.04938679789 17.1895510009+69.038358928i;8.21672845289 17.610307475+70.2348047678i;8.38754903495 18.0404968862+71.453381195i;8.56192086874 18.4806437697+72.6945717862i;8.73991778256 18.9313051512+73.9588442218i;8.92161513952 19.3930695216+75.2466477525i;9.10708986948 19.8665555936+76.5584107735i;9.29642050165 20.3524108547+77.8945385272i;9.48968719778 20.8513099354+79.2554109523i;9.68697178616 21.3639528126+80.6413806933i;9.88835779621 21.8910628682+82.0527712826i;10.0939304939 22.433384824+83.4898755046i;10.3037769178 22.9916825745+84.9529539454i;10.517985916 23.5667369381+86.4422337336i;10.7366481837 24.159343346+87.9579074706i;10.9598563014 24.7703094881+89.5001323503i;11.1877047746 25.4004529335+91.0690294634i;11.4202900733 26.0505987406+92.6646832803i;11.6577106731 26.721577072+94.2871413084i;11.9000670968 27.4142208262+95.936413913i;12.1474619572 28.1293632962+97.6124742965i;12.4 28.8678358661+99.3152586247i];
            minIdx = x < dat(end,1);
            maxIdx = x > dat(1,1);
            n(i,minIdx) = dat(end,2);
            n(i,maxIdx) = dat(1,2);
            n(i,~minIdx & ~maxIdx) = interp1(dat(:,1),dat(:,2),x(~minIdx & ~maxIdx));
        case {'Gold','Au'}
            % Data set Rakic 1998: n,k 0.2066-12.4 um
            dat = [0.2066 1.25223320824+1.98591656445i;0.21089508319 1.28860501667+1.97534779229i;0.215279458439 1.32053986347+1.96597777856i;0.219754982075 1.34888336249+1.95853718525i;0.224323549012 1.3746284289+1.95347479861i;0.228987093564 1.39882611415+1.95092054945i;0.233747590255 1.42248550103+1.95067979422i;0.238607054657 1.44647562039+1.9522622682i;0.243567544248 1.47144312879+1.95494272809i;0.248631159276 1.49775668894+1.95784420876i;0.253800043655 1.52548341112+1.96003112281i;0.259076385867 1.55439629116+1.96059922867i;0.264462419893 1.58400639322+1.95875241037i;0.269960426156 1.6136108117+1.95386083963i;0.275572732488 1.64234733089+1.94549986044i;0.281301715116 1.66924858331+1.93347278078i;0.287149799665 1.69329166131+1.91782321271i;0.293119462189 1.71344305996+1.89884348931i;0.299213230217 1.72870316606+1.87708470227i;0.305433683824 1.73815864557+1.8533704822i;0.311783456723 1.74105369497+1.82881024922i;0.318265237381 1.73689001004+1.80479849309i;0.324881770154 1.72555779493+1.78297652266i;0.331635856454 1.70748407952+1.76512629255i;0.338530355932 1.68376050988+1.75296935729i;0.345568187689 1.65618696109+1.74786631426i;0.352752331512 1.62715563874+1.7504582283i;0.360085829136 1.5993263363+1.76035121563i;0.367571785533 1.57512081975+1.77598261657i;0.375213370223 1.55616625486+1.79477454227i;0.383013818621 1.5428732988+1.8135693006i;0.390976433402 1.53428802988+1.82921617912i;0.399104585903 1.5282356392+1.83912927233i;0.407401717548 1.52166602777+1.84168456543i;0.415871341306 1.51107621084+1.83641651484i;0.42451704318 1.49291226063+1.82404632603i;0.433342483722 1.46390381788+1.80640361052i;0.442351399584 1.42132874128+1.78630032214i;0.451547605104 1.36323722341+1.76739565371i;0.460934993913 1.28868530644+1.75405729533i;0.470517540592 1.19803128159+1.75117198333i;0.480299302349 1.0933056197+1.76378471582i;0.49028442074 0.97852316585+1.79639799686i;0.50047712342 0.859575063526+1.85189192283i;0.510881725935 0.743288882063+1.93046771509i;0.521502633549 0.635790656008+2.02941584261i;0.532344343107 0.541093350059+2.14408290908i;0.543411444943 0.460739155424+2.26939106271i;0.554708624818 0.394373918234+2.40097284195i;0.56624066591 0.340593024822+2.53559991916i;0.578012450835 0.297605602816+2.67112775467i;0.590028963715 0.263621057179+2.80625688674i;0.602295292289 0.237028009053+2.94028330861i;0.614816630068 0.21645383529+3.07289808958i;0.627598278532 0.200765474753+3.204042941i;0.640645649374 0.189045341648+3.33381143007i;0.653964266795 0.180559280841+3.46238342848i;0.667559769839 0.174724390804+3.58998261813i;0.681437914784 0.171079979374+3.71684963608i;0.695604577575 0.169262728402+3.84322568932i;0.710065756315 0.168986132123+3.96934309571i;0.724827573806 0.170023837191+4.09542033789i;0.739896280135 0.172196353782+4.22165999327i;0.755278255329 0.175360582261+4.34824843716i;0.770980012047 0.179401637013+4.47535658348i;0.787008198346 0.184226512035+4.60314118084i;0.803369600489 0.189759203016+4.73174635512i;0.820071145824 0.195936968685+4.86130520764i;0.837119905712 0.202707475388+4.99194135685i;0.854523098523 0.2100266216+5.12377036316i;0.872288092693 0.217856882985+5.25690100961i;0.890422409845 0.22616605438+5.39143643137i;0.908933727968 0.234926293578+5.52747509882i;0.927829884674 0.244113394196+5.66511166545i;0.947118880514 0.253706232338+5.80443769435i;0.966808882365 0.263686345098+5.94554227814i;0.986908226887 0.274037609153+6.08851256644i;1.00742542406 0.284745995453+6.23343421379i;1.02836916076 0.295799381814+6.38039175994i;1.0497483045 0.307187409657+6.52946895243i;1.07157190709 0.318901374431+6.68074902054i;1.09384920857 0.330934141782+6.8343149079i;1.11658964103 0.343280083388+6.99024947049i;1.13980283269 0.355935027838+7.14863564526i;1.16349861191 0.368896223022+7.30955659417i;1.18768701137 0.382162307315+7.47309582753i;1.21237827235 0.395733287491+7.63933730999i;1.237582849 0.409610521811+7.80836555206i;1.26331141285 0.42379670708+7.98026568947i;1.28957485725 0.438295868801+8.15512355256i;1.31638430203 0.453113353762+8.33302572727i;1.34375109819 0.468255824603+8.51405960937i;1.37168683271 0.48373125601+8.69831345306i;1.40020333347 0.499548932349+8.88587641506i;1.42931267422 0.515719446594+9.07683859508i;1.45902717974 0.532254700497+9.27129107334i;1.48935943101 0.549167905985+9.46932594594i;1.52032227059 0.566473587824+9.67103635839i;1.55192880799 0.584187587597+9.87651653797i;1.58419242529 0.60232706909+10.085861825i;1.61712678277 0.620910525167+10.2991687034i;1.65074582469 0.639957786243+10.5165348313i;1.6850637852 0.659490030471+10.7380590705i;1.72009519439 0.679529795737+10.9638415166i;1.75585488439 0.700100993588+11.1939835286i;1.7923579957 0.721228925198+11.4285877588i;1.82961998359 0.742940299461+11.6677581822i;1.86765662461 0.765263253341+11.9116001269i;1.90648402331 0.788227374534+12.1602203034i;1.94611861905 0.81186372657+12.4137268346i;1.98657719294 0.836204876404+12.6722292853i;2.02787687497 0.861284924597+12.9358386922i;2.07003515123 0.887139538148+13.2046675927i;2.11306987137 0.913805986037+13.4788300544i;2.15699925609 0.941323177553+13.7584417033i;2.2018419049 0.969731703453+14.0436197522i;2.24761680399 0.999073880003+14.3344830281i;2.29434333425 1.02939379595+14.6311519989i;2.34204127949 1.06073736249+14.9337487993i;2.39073083481 1.09315236616+15.242397256i;2.44043261516 1.12668852494+15.557222911i;2.49116766405 1.16139754724+15.8783530446i;2.54295746248 1.19733319414+16.205916696i;2.59582393804 1.23455134469+16.5400446837i;2.64978947414 1.27311006436+16.8808696226i;2.70487691958 1.31306967671+17.2285259404i;2.76110959812 1.3544928382+17.5831498916i;2.81851131845 1.3974446162+17.9448795695i;2.8771063842 1.44199257024+18.3138549153i;2.93691960427 1.48820683645+18.6902177253i;2.99797630331 1.53616021518+19.0741116548i;3.06030233246 1.58592826192+19.4656822198i;3.1239240803 1.63758938132+19.8650767947i;3.188868484 1.69122492445+20.2724446074i;3.25516304072 1.7469192892+20.6879367303i;3.32283581931 1.80476002384+21.1117060677i;3.39191547211 1.86483793364+21.5439073394i;3.46243124716 1.92724719063+21.9846970594i;3.53441300053 1.99208544622+22.4342335112i;3.60789120897 2.05945394691+22.8926767167i;3.68289698284 2.12945765279+23.3601884018i;3.75946207926 2.20220535878+23.8369319549i;3.83761891556 2.27780981863+24.3230723814i;3.91740058299 2.35638787144+24.818776251i;3.99884086078 2.43806057054+25.3242116392i;4.08197423038 2.52295331476+25.8395480625i;4.1668358901 2.6111959817+26.3649564062i;4.25346176999 2.70292306299+26.9006088451i;4.34188854708 2.79827380114+27.4466787567i;4.43215366088 2.89739232795+28.0033406265i;4.52429532923 3.00042780404+28.5707699445i;4.61835256454 3.10753455933+29.1491430937i;4.71436519021 3.218872234+29.738637229i;4.81237385758 3.33460591977+30.3394301469i;4.91242006309 3.45490630092+30.9517001453i;5.01454616587 3.57994979477+31.5756258729i;5.11879540566 3.70991869101+32.2113861689i;5.22521192114 3.84500128958+32.8591598904i;5.3338407686 3.98539203627+33.5191257296i;5.44472794102 4.13129165578+34.191462019i;5.55792038754 4.28290728129+34.8763465244i;5.67346603336 4.44045258009+35.5739562259i;5.79141379999 4.60414787441+36.2844670861i;5.91181362602 4.77422025668+37.0080538055i;6.0347164882 4.95090369846+37.7448895647i;6.16017442306 5.13443915202+38.4951457535i;6.28824054896 5.32507464379+39.2589916857i;6.41896908852 5.52306535848+40.0365943015i;6.55241539166 5.72867371303+40.8281178548i;6.68863595894 5.94216941914+41.633723588i;6.82768846557 6.16382953324+42.4535693921i;6.96963178576 6.39393849275+43.2878094535i;7.11452601772 6.63278813733+44.1365938879i;7.26243250901 6.88067771379+45.0000683603i;7.41341388261 7.1379138633+45.8783736927i;7.56753406337 7.40481058966+46.7716454592i;7.72485830511 7.68168920691+47.6800135698i;7.88545321822 7.96887826518+48.6036018426i;8.04938679789 8.26671345301+49.5425275664i;8.21672845289 8.57553747479+50.4969010538i;8.38754903495 8.8956999019+51.4668251858i;8.56192086874 9.22755699586+52.4523949503i;8.73991778256 9.57147150229+53.4536969744i;8.92161513952 9.927812414+54.4708090531i;9.10708986948 10.296954702+55.5037996756i;9.29642050165 10.6792790132+56.5527275519i;9.48968719778 11.0751713331+57.6176411403i;9.68697178616 11.4850226129+58.6985781795i;9.88835779621 11.9092283598+59.7955652274i;10.0939304939 12.3481881894+60.9086172084i;10.3037769178 12.8023053398+62.0377369732i;10.517985916 13.2719861468+63.1829148725i;10.7366481837 13.7576394791+64.3441283502i;10.9598563014 14.2596761352+65.5213415561i;11.1877047746 14.7785082+66.7145049842i;11.4202900733 15.3145483628+67.923555139i;11.6577106731 15.8682091977+69.1484142333i;11.9000670968 16.4399024059+70.388989921i;12.1474619572 17.0300380225+71.6451750702i;12.4 17.6390235898+72.9168475779i];
            minIdx = x < dat(end,1);
            maxIdx = x > dat(1,1);
            n(i,minIdx) = dat(end,2);
            n(i,maxIdx) = dat(1,2);
            n(i,~minIdx & ~maxIdx) = interp1(dat(:,1),dat(:,2),x(~minIdx & ~maxIdx));
        case {'Silver','Ag'}
            % Data set Rakic 1998: n,k 0.2066-12.4 um
            dat = [0.2066 0.503999788208+0.870386410508i;0.21089508319 0.554820958015+0.965578874806i;0.215279458439 0.611868931581+1.0402424751i;0.219754982075 0.670439250054+1.0963412968i;0.224323549012 0.726540962221+1.13584089596i;0.228987093564 0.776528284552+1.16100205284i;0.233747590255 0.816935367231+1.17463650616i;0.238607054657 0.844330284417+1.18045934249i;0.243567544248 0.855169339166+1.18390627144i;0.248631159276 0.846331283615+1.19421880639i;0.253800043655 0.819421746645+1.22785626058i;0.259076385867 0.793290349699+1.30611642846i;0.264462419893 0.812406498312+1.43145169912i;0.269960426156 0.917874624205+1.5608011454i;0.275572732488 1.10452403555+1.62600604918i;0.281301715116 1.3192280026+1.57931089731i;0.287149799665 1.49071987712+1.41816898493i;0.293119462189 1.56281898694+1.1823642914i;0.299213230217 1.51422846209+0.932180210113i;0.305433683824 1.35789297802+0.722775908062i;0.311783456723 1.1252243269+0.589889800453i;0.318265237381 0.850147360074+0.554446352656i;0.324881770154 0.578743734444+0.636611154842i;0.331635856454 0.389430467307+0.809939685003i;0.338530355932 0.290153837234+0.988366492091i;0.345568187689 0.236246300885+1.14299721546i;0.352752331512 0.203455932597+1.27754650773i;0.360085829136 0.181911029181+1.39792320109i;0.367571785533 0.167145320914+1.50821539365i;0.375213370223 0.156810410983+1.61110819046i;0.383013818621 0.149517310499+1.70843529341i;0.390976433402 0.144378231562+1.80151065194i;0.399104585903 0.14079705276+1.89131464039i;0.407401717548 0.138360408352+1.97860222952i;0.415871341306 0.13677535281+2.06396934884i;0.42451704318 0.135831121397+2.14789577028i;0.433342483722 0.135374458563+2.23077408808i;0.442351399584 0.135293095369+2.31293005359i;0.451547605104 0.135504373517+2.39463729712i;0.460934993913 0.135947250214+2.47612826339i;0.470517540592 0.136576592738+2.55760250569i;0.480299302349 0.137359060305+2.63923308359i;0.49028442074 0.138270105886+2.72117156499i;0.50047712342 0.139291778649+2.80355198049i;0.510881725935 0.140411104466+2.88649397799i;0.521502633549 0.141618886884+2.97010535882i;0.532344343107 0.142908815745+3.05448413034i;0.543411444943 0.144276801957+3.13972017749i;0.554708624818 0.145720479208+3.22589663155i;0.56624066591 0.14723882938+3.31309099735i;0.578012450835 0.148831899984+3.40137608619i;0.590028963715 0.150500590304+3.49082079237i;0.602295292289 0.152246489084+3.58149074276i;0.614816630068 0.154071751047+3.67344884315i;0.627598278532 0.155979002845+3.7667557403i;0.640645649374 0.15797127141+3.86147021478i;0.653964266795 0.160051929504+3.95764951683i;0.667559769839 0.162224654538+4.05534965506i;0.681437914784 0.164493397713+4.15462564603i;0.695604577575 0.16686236126+4.25553173119i;0.710065756315 0.169335982084+4.35812156644i;0.724827573806 0.171918920521+4.46244838878i;0.739896280135 0.174616053216+4.5685651635i;0.755278255329 0.177432469347+4.67652471505i;0.770980012047 0.180373469619+4.7863798438i;0.787008198346 0.183444567534+4.89818343104i;0.803369600489 0.186651492603+5.01198853364i;0.820071145824 0.190000195192+5.12784847009i;0.837119905712 0.193496852774+5.2458168989i;0.854523098523 0.197147877428+5.36594789053i;0.872288092693 0.200959924413+5.48829599367i;0.890422409845 0.204939901727+5.61291629667i;0.908933727968 0.20909498055+5.73986448466i;0.927829884674 0.213432606509+5.86919689304i;0.947118880514 0.217960511702+6.00097055769i;0.966808882365 0.22268672745+6.13524326238i;0.986908226887 0.227619597745+6.27207358378i;1.00742542406 0.232767793366+6.41152093422i;1.02836916076 0.238140326666+6.55364560265i;1.0497483045 0.243746567006+6.69850879392i;1.07157190709 0.24959625685+6.84617266664i;1.09384920857 0.25569952852+6.99670036976i;1.11658964103 0.262066921616+7.15015607809i;1.13980283269 0.268709401126+7.30660502679i;1.16349861191 0.275638376222+7.46611354511i;1.18768701137 0.282865719789+7.6287490893i;1.21237827235 0.290403788672+7.79458027498i;1.237582849 0.298265444708+7.9636769089i;1.26331141285 0.306464076529+8.13611002019i;1.28957485725 0.315013622189+8.31195189123i;1.31638430203 0.323928592641+8.49127608816i;1.34375109819 0.333224096088+8.67415749101i;1.37168683271 0.34291586325+8.86067232358i;1.40020333347 0.353020273578+9.05089818309i;1.42931267422 0.363554382448+9.24491406954i;1.45902717974 0.37453594939+9.4428004149i;1.48935943101 0.385983467363+9.64463911213i;1.52032227059 0.397916193147+9.85051354396i;1.55192880799 0.41035417887+10.0605086115i;1.58419242529 0.42331830473+10.2747107629i;1.61712678277 0.436830312947+10.4932080212i;1.65074582469 0.450912842992+10.716090013i;1.6850637852 0.46558946815+10.9434479958i;1.72009519439 0.480884733448+11.1753748861i;1.75585488439 0.496824195013+11.4119652867i;1.7923579957 0.513434460902+11.6533155137i;1.82961998359 0.530743233453+11.8995236238i;1.86765662461 0.548779353227+12.1506894408i;1.90648402331 0.567572844566+12.4069145815i;1.94611861905 0.587154962849+12.6683024819i;1.98657719294 0.607558243483+12.9349584225i;2.02787687497 0.628816552695+13.2069895534i;2.07003515123 0.650965140174+13.4845049182i;2.11306987137 0.674040693628+13.7676154783i;2.15699925609 0.698081395311+14.0564341359i;2.2018419049 0.723126980572+14.3510757562i;2.24761680399 0.749218798497+14.6516571889i;2.29434333425 0.776399874689+14.9582972893i;2.34204127949 0.804714976256+15.2711169376i;2.39073083481 0.834210679055+15.5902390578i;2.44043261516 0.86493543726+15.9157886351i;2.49116766405 0.8969396553+16.2478927319i;2.54295746248 0.930275762233+16.586680503i;2.59582393804 0.964998288602+16.9322832086i;2.64978947414 1.00116394583+17.2848342262i;2.70487691958 1.03883170822+17.6444690604i;2.76110959812 1.07806289757+18.011325351i;2.81851131845 1.11892127047+18.3855428791i;2.8771063842 1.16147310837+18.7672635706i;2.93691960427 1.20578731039+19.1566314978i;2.99797630331 1.2519354889+19.5537928782i;3.06030233246 1.29999206797+19.9588960707i;3.1239240803 1.3500343847+20.3720915685i;3.188868484 1.40214279337+20.7935319897i;3.25516304072 1.45640077249+21.2233720631i;3.32283581931 1.51289503483+21.6617686123i;3.39191547211 1.57171564019+22.1088805344i;3.46243124716 1.63295611116+22.5648687752i;3.53441300053 1.69671355174+23.0298963006i;3.60789120897 1.76308876868+23.5041280626i;3.68289698284 1.83218639568+23.9877309606i;3.75946207926 1.9041150203+24.4808737984i;3.83761891556 1.97898731348+24.9837272343i;3.91740058299 2.05692016158+25.4964637269i;3.99884086078 2.13803480095+26.0192574734i;4.08197423038 2.22245695475+26.552284343i;4.1668358901 2.31031697193+27.0957218022i;4.25346176999 2.40174996825+27.6497488336i;4.34188854708 2.49689596911+28.2145458474i;4.43215366088 2.59590005391+28.7902945846i;4.52429532923 2.69891250186+29.3771780121i;4.61835256454 2.80608893878+29.9753802094i;4.71436519021 2.9175904847+30.5850862462i;4.81237385758 3.03358390194+31.2064820503i;4.91242006309 3.15424174321+31.8397542665i;5.01454616587 3.2797424994+32.485090104i;5.11879540566 3.4102707466+33.1426771747i;5.22521192114 3.54601729182+33.8127033195i;5.3338407686 3.68717931696+34.4953564235i;5.44472794102 3.83396052043+35.1908242194i;5.55792038754 3.98657125569+35.8992940792i;5.67346603336 4.14522866629+36.6209527927i;5.79141379999 4.31015681644+37.3559863338i;5.91181362602 4.48158681645+38.1045796133i;6.0347164882 4.6597569423+38.8669162182i;6.16017442306 4.84491274827+39.6431781381i;6.28824054896 5.03730717188+40.4335454766i;6.41896908852 5.23720063004+41.2381961492i;6.55241539166 5.44486110543+42.0573055671i;6.68863595894 5.66056422192+42.8910463065i;6.82768846557 5.88459330802+43.7395877636i;6.96963178576 6.11723944697+44.603095795i;7.11452601772 6.3588015123+45.481732345i;7.26243250901 6.60958618758+46.3756550575i;7.41341388261 6.86990796887+47.2850168756i;7.56753406337 7.14008914852+48.209965627i;7.72485830511 7.42045977894+49.1506435972i;7.88545321822 7.71135761469+50.1071870902i;8.04938679789 8.01312803146+51.0797259781i;8.21672845289 8.32612392049+52.0683832401i;8.38754903495 8.65070555663+53.0732744916i;8.56192086874 8.98724043875+54.0945075055i;8.73991778256 9.33610310082+55.1321817261i;8.92161513952 9.69767489212+56.1863877775i;9.10708986948 10.0723437251+57.2572069682i;9.29642050165 10.4605037896+58.3447107935i;9.48968719778 10.8625552319+59.4489604375i;9.68697178616 11.2789037969+60.5700062781i;9.88835779621 11.7099604336+61.7078873953i;10.0939304939 12.1561408611+62.8626310875i;10.3037769178 12.6178650954+64.0342523968i;10.517985916 13.0955569356+65.222753648i;10.7366481837 13.5896434094+66.428124002i;10.9598563014 14.1005541777+67.6503390297i;11.1877047746 14.6287208969+68.8893603073i;11.4202900733 15.1745765412+70.1451350383i;11.6577106731 15.738554683+71.4175957048i;11.9000670968 16.3210887347+72.7066597524i;12.1474619572 16.92261115+74.012229313i;12.4 17.5435525892+75.3341909676i];
            minIdx = x < dat(end,1);
            maxIdx = x > dat(1,1);
            n(i,minIdx) = dat(end,2);
            n(i,maxIdx) = dat(1,2);
            n(i,~minIdx & ~maxIdx) = interp1(dat(:,1),dat(:,2),x(~minIdx & ~maxIdx));
        case 'none'
            % This will be used for pinholes
            n(i,:) = 1;
        otherwise
            error('RefactiveIndex:unknownMaterial', 'Unknown material')
    end
end

end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%