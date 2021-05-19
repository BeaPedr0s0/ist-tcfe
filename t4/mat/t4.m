%gain stage

VT=25e-3;
BFN=178.7;
VAFN=69.7;
RE1=100;
RC1=1000;
RB1=80000;
RB2=20000;
VBEON=0.7;
VCC=12;
RS=100;

fich = fopen("initial_values_tab.tex","w")
	fprintf(fich, "VT & %f \\\\ \\hline \n", VT)
	fprintf(fich, "BFN& %f \\\\ \\hline \n", BFN)
	fprintf(fich, "VAFN& %f \\\\ \\hline \n", VAFN)
	fprintf(fich, "RE1& %f \\\\ \\hline \n", RE1)
	fprintf(fich, "RC1& %f \\\\ \\hline \n", RC1)
	fprintf(fich, "RB1& %f \\\\ \\hline \n", RB1)
	fprintf(fich, "RB2& %f \\\\ \\hline \n", RB2)
	fprintf(fich, "VBEON& %f \\\\ \\hline \n", VBEON)
	fprintf(fich, "VCC& %f \\\\ \\hline \n", VCC)
	fprintf(fich, "RS& %f \\\\ \\hline \n", RS)
fclose(fich)


RB=1/(1/RB1+1/RB2);
VEQ=RB2/(RB1+RB2)*VCC;
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1);
IC1=BFN*IB1;
IE1=(1+BFN)*IB1;
VE1=RE1*IE1;
VO1=VCC-RC1*IC1;
VCE=VO1-VE1;


gm1=IC1/VT;
rpi1=BFN/gm1;
ro1=VAFN/IC1;

AV1 = RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2);

AV1simple = gm1*RC1/(1+gm1*RE1);

RE1=0;
AV1 = RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2);
AV1simple = gm1*RC1/(1+gm1*RE1);

RE1=100;

ZI1 = ((ro1+RC1+RE1)*(RB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1);

ZX = ro1*((RB+rpi1)*RE1/(RB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RB)+1/RE1+gm1*rpi1/(rpi1+RB)));

ZO1 = 1/(1/ZX+1/RC1);

tab=fopen("gain.tex","w");
fprintf(tab, "Gain & %f \\\\ \\hline \n", AV1simple);
fprintf(tab, "Input Impedance & %fk \\\\ \\hline \n", ZI1/1000);
fprintf(tab, "Output Impedance & %f \\\\ \\hline \n", ZO1);
fclose(tab);

%ouput stage
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2)

ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)

ZO2 = 1/(gm2+gpi2+go2+ge2)


tab=fopen("output.tex","w");
fprintf(tab, "Gain & %f \\\\ \\hline \n", AV2);
fprintf(tab, "Input Impedance & %fk \\\\ \\hline \n", ZI2/1000);
fprintf(tab, "Output Impedance & %f \\\\ \\hline \n", ZO2);
fclose(tab);