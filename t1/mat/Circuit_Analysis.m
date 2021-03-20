
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

printf ("malhas_TAB\n");
printf ("Ia = %e \n", C(1));
printf ("Ib = %e \n", C(2));
printf ("Ic = %e \n", C(3));
printf ("malhas_END\n");



printf ("nos_TAB\n");
printf ("V0 = %e \n", F(1));
printf ("V1 = %e \n", F(2));
printf ("V2 = %e \n", F(3));
printf ("V3 = %e \n", F(4));
printf ("V4 = %e \n", F(5));
printf ("V5 = %e \n", F(6));
printf ("V6 = %e \n", F(7));
printf ("V7 = %e \n", F(8));
printf ( "nos_END\n");


Ia=(F(2)-F(1))/R1;
Ib=(F(3)-F(2))/R2;
Ic=(F(4)-F(7))/R6;

printf ("confirmacao_TAB\n");
printf ("Ia = %e \n", Ia);
printf ("Ib = %e \n", Ib);
printf ("Ic = %e \n", Ic);
printf ("confirmacao_END\n");

