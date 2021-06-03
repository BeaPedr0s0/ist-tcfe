%Lab4

C1 =500e-06;
C2 =500e-06;
RB1 =20e3 ;
RB2 = 2e3;
RC1 = 3e3;
RE1 =0.1e3;
Rout=0.2e3;
Cout=200e-6;
RL=8;

%Stage1
VT=25e-3;
BFN=178.7;
VAFN=69.7;
VBEON=0.7;
VCC=12;
Rin=100;

%Stage2
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7

printf ("valores_intro_TAB\n");
printf ("Cin = %e \n", C1);
printf ("CE = %e \n", C2);
printf ("Cout = %e \n", Cout);
printf ("R1 = %e \n", RB1);
printf ("R2 = %e \n", RB2);
printf ("RC = %e \n", RC1);
printf ("RE = %e \n", RE1);
printf ("Rout = %e \n", RE1);
printf ("Vcc = %e \n", VCC);
printf ("valores_intro_END\n\n");

%Values for ngspice
%FILE1

file1=fopen("values1.cir",'w');
fprintf(file1, ".OP\n\n");
fprintf(file1, "Vcc vcc 0 12\n");
fprintf(file1, "Vin in 0 0 \n");
fprintf(file1, "Rin in in2 100 \n\n");
fprintf(file1,"*input  coupling capacitor\n");
fprintf(file1, "Ci in2 base 500u \n\n");
fprintf(file1,"*bias circuit\n");
fprintf(file1, "R1 vcc base 20k \n");
fprintf(file1, "R2 base 0 2k \n\n");
fprintf(file1,"*gain stage\n");
fprintf(file1, "Q1 coll base emit BC547A\n");
fprintf(file1, "Rc vcc coll 3k\n");
fprintf(file1, "Re emit 0 0.1k\n\n");
fprintf(file1,"*bypass capacitor\n");
fprintf(file1, "Cb emit 0 500u\n\n");
fprintf(file1,"*output stage\n");
fprintf(file1, "Q2 0 coll emit2 BC557A\n");
fprintf(file1, "Rout emit2 vcc 0.2k\n\n");
fprintf(file1,"*output coupling capacitor\n");
fprintf(file1, "Co emit2 out 200u\n\n");
fprintf(file1,"*fonte de teste\n");
fprintf(file1, "VL out 0 ac 1.0 sin(0 10m 1k)\n\n");
fprintf(file1, ".END\n\n");
fclose(file1);


%%%FILE2
file2=fopen("description.cir",'w');
fprintf(file2, ".OP\n\n");
fprintf(file2, "Vcc vcc 0 12 \n");
fprintf(file2, "Vin in 0 0 ac 1.0 sin(0 10m 1k) \n");
fprintf(file2, "Rin in in2 100\n\n");
fprintf(file2,"*input  coupling capacitor\n");
fprintf(file2, "Ci in2 base 500u\n\n");
fprintf(file2,"*bias circuit\n");
fprintf(file2, "R1 vcc base 20k \n");
fprintf(file2, "R2 base 0 2k \n\n");
fprintf(file2,"*gain stage\n");
fprintf(file2, "Q1 coll base emit BC547A\n");
fprintf(file2, "Rc vcc coll 3k\n");
fprintf(file2, "Re emit 0 0.1k\n\n");
fprintf(file2,"*bypass capacitor\n");
fprintf(file2, "Cb emit 0 500u\n\n");
fprintf(file2,"*output stage\n");
fprintf(file2, "Q2 0 coll emit2 BC557A\n");
fprintf(file2, "Rout emit2 vcc 0.2k\n\n");
fprintf(file2,"*output coupling capacitor\n");
fprintf(file2, "Cout emit2 out 200u\n\n");
fprintf(file2,"*load\n");
fprintf(file2, "RL out 0 8\n\n");
fprintf(file2, ".END\n\n");
fclose(file2);


%//////////////////GAIN STAGE/////////////////////%
RB=1/(1/RB1+1/RB2);
VEQ=RB2/(RB1+RB2)*VCC;
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1);
IC1=BFN*IB1;
IE1=(1+BFN)*IB1;
VE1=RE1*IE1;
VO1=VCC-RC1*IC1;
VCE=VO1-VE1;

gm1 = IC1/VT;
rpi1 = BFN/gm1;
ro1 = VAFN/IC1;

RinB=RB*Rin/(RB+Rin);

AV1 = RinB/Rin * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RinB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2);
AV1db = 20*log10(abs(AV1));

Zin1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)));
Zoutut1 = 1/(1/ro1+1/RC1);

%//////////////////OUTPUT STAGE/////////////////////%
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2);
AV2db= 20*log10(abs(AV2));

Zin2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);
Zoutut2 = 1/(gm2+gpi2+go2+ge2);

%total
gB = 1/(1/gpi2+Zoutut1);
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1;
AV_TOT_db = 20*log10(abs(AV));

ZI=Zin1;
Zout=1/(go2+gm2/gpi2*gB+ge2+gB);

%LowCutOff Frequency
R1S = Rin + (1/(1/RB + 1/rpi1));
R2S = RL + (1/(1/RC1 + 1/ro1));
R3S = 1/((1/RE1) + (1/(rpi1 + (1/(1/Rin + 1/RB)))) + ((gm1*rpi1)/(rpi1 + (1/(1/Rin + 1/RB)))));
wL = 1/(R1S*C1) + 1/(R2S*C2) + 1/(R3S*Cout);
LowFreq = wL/(2*pi);

%HighCutOff Frequency
Cpi = 16.1e-12;
Co = 4.388e-12;
wH = 1/(Cpi*rpi1 + Co*ro1);
HighFreq = wH/(2*pi);

w = logspace(1,12);
Tdb = ones(1,length(w));
for k = 1:length(w)
	T = 10^((AV_TOT_db-20*log10(wL))/20)*(w(k)/(w(k)/wL + 1))*(1/(w(k)/wH + 1));
	Tdb(k) = 20*log10(abs(T));
end

cost = 1e-3*(RE1 + RC1 + RB1 + RB2 + RE2) + 1e6*(C1 + C2 + Cout) + 2*0.1;
Merit = (abs(AV_TOT_db) * (HighFreq-LowFreq))/(cost * LowFreq);

printf ("cost= %e \n", cost);

%%Final Result Tables and Plots
%%%%%%PONTO1%%%%%%%
printf ("ponto1_TAB\n");
printf ("IB1 = %e \n", IB1);
printf ("IC1 = %e \n", IC1);
printf ("IE1 = %e \n", IE1);
printf ("VColl = %e \n", VO1);
printf ("VBase = %e \n", VEQ);
printf ("VEmit = %e \n", VE1);
printf ("VCE = %e V\n", VCE);
printf ("VBEON = %e V \n", VBEON);
printf ("VEC = %e V\n", VO2);
printf ("VEBON = %e V \n", VEBON);
printf ("IB2 = %e A \n", IC2-IE2);
printf ("IC2 = %e A \n", IC2);
printf ("IE2 = %e A \n", IE2);
printf ("ponto1_END\n\n");

%%%%%%PONTO2%%%%%%%
printf ("Z_TAB\n");
printf ("Zin1 = %e Ohm \n", real(Zin1));
printf ("Zoutut1 = %e Ohm \n", Zoutut1);
printf ("Zin2 = %e Ohm \n", Zin2);
printf ("Zoutut2 = %e Ohm \n", Zoutut2);
printf ("Zout = %e \n", Zout);
printf ("Z_END\n\n");

printf ("r_theo_TAB\n");
printf ("Total Gain (dB)  = %e V\n", AV_TOT_db);
printf ("Low Cut Off Frequency= %e Hz \n", LowFreq);
printf ("High Cut Off Frequency= %e Hz \n", HighFreq);
printf ("Bandwidth= %e Hz \n", HighFreq-LowFreq);
printf ("Cost = %e MU's\n", cost);
printf ("Merit = %e \n", Merit);
printf ("r_theo_END\n\n");


%%%%%%Ponto3%%%%%%%
theo = figure ();
plot(log10(w/(2*pi)),Tdb,"g");
legend("v_o(f)/v_i(f)");
xlabel ("Log10(Frequency [Hz])");
ylabel ("Gain");
print (theo, "theo", "-depsc");




