
R1 = 1.04111259479E3;
R2 = 2.09945227782E3;
R3 = 3.13109125645E3;
R4 = 4.11947040212E3;
R5 = 3.1155879392E3;
R6 = 2.04799381798E3; 
R7 = 1.02754401839E3;

G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

Va = 5.06871572779;
Id = 1.04127523824E-3;
Kb = 7.28747116393E-3;
Kc = 8.11568444746E-3;

A=[R1+R4+R3, -R3, -R4; -R4, 0, R6+R7-Kc+R4;-R3*Kb, R3*Kb-1, 0];
B=[-Va; 0; 0];

C=inv(A)*B;

Ia= C(1);
Ib= C(2);
Ic= C(3);

D =[1 0 0 -1 0 0 0 0;
    -G1 G1+G2+G3 -G2 0 -G3 0 0 0;
    0 -Kb 0 0 Kb+G5 -G5 0 0;
    0 0 0 G6 0 0 -G6-G7 G7;
    1 0 0 0 0 0 0 0;
    0 Kb+G2 -G2 0 -Kb 0 0 0;
    0 G1 0 -G4-G6 G4 0 G6 0;
    0 0 0 -Kc*G6 1 0 Kc*G6 -1];

E=[Va; 0 ; -Id; 0 ; 0 ; 0 ; 0; 0];
F=inv(D)*E;
V0=F(1);
V1=F(2);
V2=F(3);
V3=F(4);
V4=F(5);
V5=F(6);
V6=F(7);
V7=F(8);


file=fopen("tabelamalhas.tex",'w');

fprintf(file,'%s', 'ola');
%fprintf(file,"Método das Malhas: \n");
%fprintf(file,"Ia=%f \n",Ia);
%fprintf(file,"Ib=%f \n",Ib);
%fprintf(file,"Ic=%f \n",Ic);
%fprintf(file,"Id=%f \n\n",Id);

fclose(file);

file1=fopen("tabelanos.tex",'w');

fprintf(file,"Método dos Nós: \n");
fprintf(file,"V0=%f \n",V0);
fprintf(file,"V1=%f \n",V1);
fprintf(file,"V2=%f \n",V2);
fprintf(file,"V3=%f \n",V3);
fprintf(file,"V4=%f \n",V4);
fprintf(file,"V5=%f \n",V5);
fprintf(file,"V6=%f \n",V6);
fprintf(file,"V7=%f",V7);

fclose(file1);