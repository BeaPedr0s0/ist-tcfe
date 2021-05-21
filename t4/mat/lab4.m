%//////////////////GAIN STAGE/////////////////////%

%%%%%%% optimized variables %%%%%%%%%
C1 = 1.9e-06;
C2 = 1.0e-06;
RB1 = 76000;
RB2 = 19000;
RC1 = 20000;
RE1 = 3900;
%%%%%%%%

i = complex(0,1);
f = logspace(1,8, 100);
w = 2*pi*f;

vi = 0.01;
VT=25e-3;
BFN=178.7;
VAFN=69.7;
VBEON=0.7;
VCC=12;
RS=100;

%Values for ngspice
%FILE1
file1=fopen("values1.cir",'w');
fprintf(file1, ".OP\n\n");
fprintf(file1, "Vcc vcc 0 12\n");
fprintf(file1, "Vin in 0 0 \n");
fprintf(file1, "Rin in in2 100 \n\n");
fprintf(file1,"*input  coupling capacitor\n");
fprintf(file1, "Ci in2 base 1.9345e-06 \n\n");
fprintf(file1,"*bias circuit\n");
fprintf(file1, "R1 vcc base 76200 \n");
fprintf(file1, "R2 base 0 19401 \n\n");
fprintf(file1,"*gain stage\n");
fprintf(file1, "Q1 coll base emit BC547A\n");
fprintf(file1, "Rc vcc coll 20000\n");
fprintf(file1, "Re emit 0 3892.8\n\n");
fprintf(file1,"*bypass capacitor\n");
fprintf(file1, "Cb emit 0 1.0345e-06\n\n");
fprintf(file1,"*output stage\n");
fprintf(file1, "Q2 0 coll emit2 BC557A\n");
fprintf(file1, "Rout emit2 vcc 100\n\n");
fprintf(file1,"*output coupling capacitor\n");
fprintf(file1, "Co emit2 out 1u\n\n");
fprintf(file1,"*fonte de teste\n");
fprintf(file1, "VL out 0 ac 1.0 sin(0 10m 1k)\n\n");
fprintf(file1, ".END\n\n");
fclose(file1);

%%%FILE2
file2=fopen("description.cir",'w');
fprintf(file1, ".OP\n\n");
fprintf(file2, "Vcc vcc 0 12 \n");
fprintf(file2, "Vin in 0 0 ac 1.0 sin(0 10m 1k) \n");
fprintf(file2, "Rin in in2 {R_init} \n\n");
fprintf(file2,"*input  coupling capacitor\n");
fprintf(file2, "Ci in2 base {C_i} \n\n");
fprintf(file2,"*bias circuit\n");
fprintf(file2, "R1 vcc base {R_1} \n");
fprintf(file2, "R2 base 0 {R_2} \n\n");
fprintf(file2,"*gain stage\n");
fprintf(file2, "Q1 coll base emit BC547A\n");
fprintf(file2, "Rc vcc coll {R_c}\n");
fprintf(file2, "Re emit 0 {R_e}\n\n");
fprintf(file2,"*bypass capacitor\n");
fprintf(file2, "Cb emit 0 {C_b}\n\n");
fprintf(file2,"*output stage\n");
fprintf(file2, "Q2 0 coll emit2 BC557A\n");
fprintf(file2, "Rout emit2 vcc {R_out}\n\n");
fprintf(file2,"*output coupling capacitor\n");
fprintf(file2, "Co emit2 out {C_out}\n\n");
fprintf(file2,"*load\n");
fprintf(file2, "RL out 0 {R_L}\n\n");
fprintf(file2, ".END\n\n");
fclose(file2);

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
Zc1 = 1./(i*w*C1);
Zc2 = 1./(i*w*C2);

I=ones(length(w),4);
voE=zeros(1,length(w));
AV1=zeros(1,length(w));

for k=1:length(w)
A = [0, -RB, 0, 0, RS+RB;
    -RE1, -Zc2(k), 0, RE1 + Zc2(k), 0; 
    RE1 + ro1 + RC1, 0, -ro1, -RE1, 0;
    0, gm1*rpi1, 1, 0, 0;
    0, RB+Zc1(k)+rpi1+Zc2(k), 0, -Zc2(k), -RB];

N = [vi; 0; 0; 0; 0];
res = A\N;

voE(k) = abs(RC1 * res(1));
AV1(k) = voE(k)/vi;
endfor

ZI1 = 1/(1/RB + 1/rpi1);
ZO1 = 1/(1/ro1 + 1/RC1);

%???????????????????????????????????
%ZI1 = ((ro1+RC1+RE1)*(RB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1);
%ZX = ro1*((RB+rpi1)*RE1/(RB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RB)+1/RE1+gm1*rpi1/(rpi1+RB)));
%ZO1 = 1/(1/ZX+1/RC1);


%//////////////////OUTPUT STAGE/////////////////////%
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

ro2 = 1/go2;
rpi2 = 1/gpi2;

I2=zeros(length(w),3);
vo2=zeros(1,length(w));
AV2 = zeros(1,length(w));

for k=1:length(w)
A2 = [rpi2+RE2, -RE2, 0;
    gm2*rpi2, 0, 1; 
    -RE2, RE2+ro2, -ro2];

N2 = [voE(k); 0; 0];
res = A2\N2;

vo2(k) = abs((res(1)-res(2))*RE2);
AV2(k) = vo2(k)/voE(k);
endfor


ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);
ZO2 = 1/(gm2+gpi2+go2+ge2);
ZO = 1/(gm2*(rpi2/(rpi2+ZO1))+(1/(rpi2+ZO1))+go2+ge2);
AV = AV1.*AV2;


%%%%%%PLOT%%%%%%%
theo = figure ();
plot (log10(f), 20*log10(AV), "g");
legend("v_o(f)/v_i(f)");
xlabel ("Frequency [Hz]");
ylabel ("Gain");
print (theo, "theo", "-depsc");
%%%%%%% end plot %%%%%%%%

AVdb = 20*log10(AV);
maximo = max(AVdb);
k = 1;

while  0.05 < ((maximo - AVdb(k))/maximo)
    k = k + 1;
endwhile

lowerCutoff = (w(k))/(2*pi);
highCutoff = 10^7;

bandwidth = highCutoff - lowerCutoff;
cost = 1e-3*(RE1 + RC1 + RB1 + RB2 + RE2) + 1e6*(C1 + C2) + 2*0.1;

AV=abs(AV);
Merit = (max(AV) * bandwidth)/(cost * lowerCutoff)

printf ("ponto2_TAB\n");
printf ("AV1 = %e \n", max(AV1));
printf ("ZI1 = %e \n", ZI1);
printf ("ZO1 = %e \n", ZO1);
printf ("AV2 = %e V \n", max(AV2));
printf ("ZI2 = %e \n", ZI2);
printf ("ZO2 = %e \n", ZO2);
printf ("ZO = %e \n", ZO);
printf ("ponto2_END\n\n");
