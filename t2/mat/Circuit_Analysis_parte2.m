format long e

printf ("\n\n-----------Parte 1-----------\n\n");

file=fopen("../data.txt",'r');
t=importdata("../data.txt",'\t',14);
R1=sscanf(char (t(4,1)), "Values: R1 = %e")*1000;
R2=sscanf(char (t(5,1)), "R2 = %e")*1000;
R3=sscanf(char (t(6,1)), "R3 = %e")*1000;
R4=sscanf(char (t(7,1)), "R4 = %e")*1000;
R5=sscanf(char (t(8,1)), "R5 = %e")*1000;
R6=sscanf(char (t(9,1)), "R6 = %e")*1000;
R7=sscanf(char (t(10,1)), "R7 = %e")*1000;
Vs=sscanf(char (t(11,1)), "Vs = %e");
C =sscanf(char (t(12,1)), "C = %e")/1000000;
Kb=sscanf(char (t(13,1)), "Kb = %e")/1000;
Kd=sscanf(char (t(14,1)), "Kd = %e\n")*1000;
fclose (file);

printf ("initial_data_TAB\n");
printf ("R1 = %e \n", R1);
printf ("R2 = %e \n", R2);
printf ("R3 = %e \n", R3);
printf ("R4 = %e \n", R4);
printf ("R5 = %e \n", R5);
printf ("R6 = %e \n", R6);
printf ("R7 = %e \n", R7);
printf ("Vs = %e \n", Vs);
printf ("C = %e \n", C);
printf ("Kb = %e \n", Kb);
printf ("Kd = %e \n", Kd);
printf ("initial_data_END\n");

% writing the data in first.cir for the first simulation in ngspice
file1=fopen("first.cir",'w');
fprintf(file1, ".OP\n\n");
fprintf(file1,"*INDEPENDENT VOLTAGE SOURCE, CONNECTED TO GROUND. AND NODE X. VALUE IN VOLTS.\n");
fprintf(file1, "Vs 1 0 %e \n\n",Vs);
fprintf(file1,"*RESISTORS, CONNECTED BETWEEN THE NODES X X. RESISTANCE VALUE IN OHM.\n");
fprintf(file1, "R1 1 2 %e \n",R1);
fprintf(file1, "R2 3 2 %e \n",R2);
fprintf(file1, "R3 2 5 %e \n",R3);
fprintf(file1, "R4 5 0 %e \n",R4);
fprintf(file1, "R5 5 6 %e \n",R5);
fprintf(file1, "R6 9 0 %e \n",R6);
fprintf(file1, "R7 7 8 %e \n\n",R7);
fprintf(file1,"*DEPENDENT CURRENT SOURCE. CONNECTED TO NODES X & X. DEPENDS OF THE VOLTAGE BETWEEN THE NODES 1 AND 4. VALUE AMPLIFIED BY CONSTANT Z. VALUE IN AMPERES.\n");
fprintf(file1, "Gb 6 3 (2,5) %e \n\n",Kb);
fprintf(file1,"*DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X . DEPENDS OF THE CURRENT THATS FLOWS IN THE VAUX SOURCE. IS AMPLIFIED BY CONSTANTE Y. VALUE IN VOLTS.\n");
fprintf(file1, "Hd 5 8 Vaux %e \n\n",Kd);
fprintf(file1,"*INDEPENDENT VOLTAGE SOURCE TO USE IN THE DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file1, "Vaux 9 7 0 \n\n");
fprintf(file1,"*CAPACITOR. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file1, "C1 6 8 %e \n\n",C);
fprintf(file1, ".END\n\n");
fclose(file1);

% Calculation of the inverse of the resistors. Useful to the Node Method
G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

% Insertion of the matrices that will allow us to calculate the voltages of 
% each node by the Node Method
A=[1 0 0 0 0 0 0 0;
    -G1 G1+G2+G3 -G2 0 -G3 0 0 0;
    0 -Kb-G2 G2 0 Kb 0 0 0;
    0 0 0 1 0 0 0 0;
    0 -G3 0 -G4 G3+G4+G5 -G5 -G7 G7;
    0 Kb 0 0 -Kb-G5 G5 0 0;
    0 0 0 -G6 0 0 G6+G7 -G7;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1];

B=[Vs; 0 ; 0; 0 ; 0 ; 0 ; 0; 0];

% A*B=C <=>
V=inv(A)*B;

Vx=V(6)-V(8);

%Currents Computation
Ia=(V(2)-V(1))/R1;
Ib=Kb*(V(2)-V(5));
Ic=(V(7)-V(8))/R7;
Id=Ib-((V(5)-V(6))/R5);


% Creating a table in the mat folder with the voltages results of the Node Method
printf ("nos_part1_TAB\n");
printf ("c1 = %e \n", 0);
printf ("Ib = %e \n", Ib);
printf ("I1 = %e \n", -Ia);
printf ("I2 = %e \n", Ib);
printf ("I3 = %e \n", Ib-Ia);
printf ("I4 = %e \n", -Ia+Ic);
printf ("I5 = %e \n", Ib-Id);
printf ("I6 = %e \n", -Ic);
printf ("I7 = %e \n", Ic);
printf ("V1 = %e \n", V(1));
printf ("V2 = %e \n", V(2));
printf ("V3 = %e \n", V(3));
printf ("V4 = %e \n", V(4));
printf ("V5 = %e \n", V(5));
printf ("V6 = %e \n", V(6));
printf ("V7 = %e \n", V(7));
printf ("V8 = %e \n", V(8));
printf ("nos_part1_END\n");

% Calculation of the inverse of the resistors. Useful to the Node Method
G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf ("\n\n-----------Parte 2-----------\n\n");

% Insertion of the matrices that will allow us to calculate the voltages of 
% each node
D =[1 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0;
    -G1 G1+G2+G3 -G2 0 -G3 0 0 0 0;
    0 -Kb-G2 G2 0 Kb 0 0 0 0;   
    0 -G1 0 0 -G4 0 -G6 0 0;
    0 0 0 0 0 0 G6+G7 -G7 0;
    0 Kb 0 0 -Kb-G5 G5 0 0 1;
    0 0 0 0 0 1 0 -1 0;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1 0];
    
E = [0; 0; 0; 0; 0; 0; 0; Vx; 0];

% D*E=F <=>
F=inverse(D)*E;

%computation of the equivalent resistor (Req=Vx/(-Ix))
Req=(Vx/(-F(9)));

%time constant (tau=Req*C)
tau=(Req*C);

% Creating a table in the mat folder with the voltages results
printf ("voltages_p2_TAB\n");
printf ("V1 = %e \n", F(1));
printf ("V2 = %e \n", F(2));
printf ("V3 = %e \n", F(3));
printf ("V4 = %e \n", F(4));
printf ("V5 = %e \n", F(5));
printf ("V6 = %e \n", F(6));
printf ("V7 = %e \n", F(7));
printf ("V8 = %e \n", F(8));
printf ("Ix = %e \n", F(9));
printf ("Req = %e \n", Req);
printf ("tau = %e \n", tau);
printf ( "voltages_p2_END\n\n");


% writing the data in second.cir for the third simulation in ngspice
file2=fopen("second.cir",'w');
fprintf(file2, ".OP\n\n");
fprintf(file2,"*INDEPENDENT VOLTAGE SOURCE, CONNECTED TO GROUND. AND NODE X. VALUE IN VOLTS.\n");
fprintf(file2, "Vs 1 0 0 \n\n");
fprintf(file2,"*RESISTORS, CONNECTED BETWEEN THE NODES X X. RESISTANCE VALUE IN OHM.\n");
fprintf(file2, "R1 2 1 %e \n",R1);
fprintf(file2, "R2 3 2 %e \n",R2);
fprintf(file2, "R3 5 2 %e \n",R3);
fprintf(file2, "R4 0 5 %e \n",R4);
fprintf(file2, "R5 5 6 %e \n",R5);
fprintf(file2, "R6 0 9 %e \n",R6);
fprintf(file2, "R7 7 8 %e \n\n",R7);
fprintf(file2,"*DEPENDENT CURRENT SOURCE. CONNECTED TO NODES X & X. DEPENDS OF THE VOLTAGE BETWEEN THE NODES 1 AND 4. VALUE AMPLIFIED BY CONSTANT Z. VALUE IN AMPERES.\n");
fprintf(file2, "Gb 6 3 (2,5) %e \n\n",Kb);
fprintf(file2,"*DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X . DEPENDS OF THE CURRENT THATS FLOWS IN THE VAUX SOURCE. IS AMPLIFIED BY CONSTANTE Y. VALUE IN VOLTS.\n");
fprintf(file2, "Hd 5 8 Vaux %e \n\n",Kd);
fprintf(file2,"*INDEPENDENT VOLTAGE SOURCE TO USE IN THE DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file2, "Vaux 9 7 0 \n\n");
fprintf(file2,"*CAPACITOR. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file2, "Vx 6 8 %e \n\n",Vx);
fprintf(file2, ".END\n\n");
fclose(file2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf ("\n\n-----------Parte 3-----------\n\n");


%time axis: 0 to 20ms with 1us steps
t = 0:(20e-3)/1000:20e-3; %s

V6_n=(F(6))*e.^(-t/tau);

%Plot natural solution 
natural = figure ();
plot (t*1000, V6_n, "r", "linewidth",4);
xlabel ("t [ms]");
ylabel ("V [V]");
h=legend({"V_{natural}"});
legend(h,"location", "northeast");
set (gca,"fontsize",16,"linewidth",2);
set (h,"fontsize",16);
print (natural, "natural.eps", "-depsc");

% writing the data in third.cir for the third simulation in ngspice
file3=fopen("third.cir",'w');
fprintf(file3, ".OP\n\n");
fprintf(file3,"*INDEPENDENT VOLTAGE SOURCE, CONNECTED TO GROUND. AND NODE X. VALUE IN VOLTS.\n");
fprintf(file3, "Vs 1 0 0 \n\n");
fprintf(file3,"*RESISTORS, CONNECTED BETWEEN THE NODES X X. RESISTANCE VALUE IN OHM.\n");
fprintf(file3, "R1 1 2 %e \n",R1);
fprintf(file3, "R2 3 2 %e \n",R2);
fprintf(file3, "R3 2 5 %e \n",R3);
fprintf(file3, "R4 5 0 %e \n",R4);
fprintf(file3, "R5 5 6 %e \n",R5);
fprintf(file3, "R6 9 0 %e \n",R6);
fprintf(file3, "R7 7 8 %e \n\n",R7);
fprintf(file3,"*DEPENDENT CURRENT SOURCE. CONNECTED TO NODES X & X. DEPENDS OF THE VOLTAGE BETWEEN THE NODES 1 AND 4. VALUE AMPLIFIED BY CONSTANT Z. VALUE IN AMPERES.\n");
fprintf(file3, "Gb 6 3 (2,5) %e \n\n",Kb);
fprintf(file3,"*DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X . DEPENDS OF THE CURRENT THATS FLOWS IN THE VAUX SOURCE. IS AMPLIFIED BY CONSTANTE Y. VALUE IN VOLTS.\n");
fprintf(file3, "Hd 5 8 Vaux %e \n\n",Kd);
fprintf(file3,"*INDEPENDENT VOLTAGE SOURCE TO USE IN THE DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file3, "Vaux 9 7 0 \n\n");
fprintf(file3,"*CAPACITOR. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file3, "C1 6 8 %e \n\n",C);
fprintf(file3,".IC v(6)=%e v(8)=%e\n\n",V(6), V(8));
fprintf(file3, ".END\n\n");
fclose(file3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf ("\n\n-----------Parte 4-----------\n\n");

%frequency [Hz]
f=1000;
w=2*pi*f;
Zc=1/(j*w*C);

G =[G1 -G1-G2-G3 G2 0 G3 0 0 0;
    0 Kb+G2 -G2 0 -Kb 0 0 0;
    0 Kb 0 0 -Kb-G5 G5+1/Zc 0 -1/Zc;
    0 0 0 -G6 0 0 G6+G7 -G7;
    0 G3 0 G4 -G3-G5-G4 G5+1/Zc G7 -G7-1/Zc;
    1 0 0 -1 0 0 0 0;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1;
    0 0 0 1 0 0 0 0];
    
H =[0; 0; 0; 0; 0; 1; 0; 0];

% G*I=H <=>
I=G\H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for a=1:length(I)
  phi(a)= arg(I(a));
endfor


%imprimir amplitude complexa
% Creating a table in the mat folder with the voltages results
printf ("voltages_p4_TAB\n");
printf ("V1 = %e \n", abs(I(1)));
printf ("V2 = %e \n", abs(I(2)));
printf ("V3 = %e \n", abs(I(3)));
printf ("V4 = %e \n", abs(I(4)));
printf ("V5 = %e \n", abs(I(5)));
printf ("V6 = %e \n", abs(I(6)));
printf ("V7 = %e \n", abs(I(7)));
printf ("V8 = %e \n", abs(I(8)));
printf ( "voltages_p4_END\n");

printf ("phase_p4_TAB\n");
printf ("Phase1 = %e \n", phi(1));
printf ("Phase2 = %e \n", phi(2));
printf ("Phase3 = %e \n", phi(3));
printf ("Phase4 = %e \n", phi(4));
printf ("Phase5 = %e \n", phi(5));
printf ("Phase6 = %e \n", phi(6));
printf ("Phase7 = %e \n", phi(7));
printf ("Phase8 = %e \n", phi(8));
printf ( "phase_p4_END\n");


% writing the data in fourth.cir for the fourth simulation in ngspice
file4=fopen("fourth.cir",'w');
fprintf(file4, ".OP\n\n");
fprintf(file4,"*INDEPENDENT VOLTAGE SOURCE, CONNECTED TO GROUND. AND NODE X. VALUE IN VOLTS.\n");
fprintf(file4, "Vs 1 0 0.0 ac 1.0 sin(0 1 1000) \n\n");
fprintf(file4,"*RESISTORS, CONNECTED BETWEEN THE NODES X X. RESISTANCE VALUE IN OHM.\n");
fprintf(file4, "R1 1 2 %e \n",R1);
fprintf(file4, "R2 3 2 %e \n",R2);
fprintf(file4, "R3 2 5 %e \n",R3);
fprintf(file4, "R4 5 0 %e \n",R4);
fprintf(file4, "R5 5 6 %e \n",R5);
fprintf(file4, "R6 9 0 %e \n",R6);
fprintf(file4, "R7 7 8 %e \n\n",R7);
fprintf(file4,"*DEPENDENT CURRENT SOURCE. CONNECTED TO NODES X & X. DEPENDS OF THE VOLTAGE BETWEEN THE NODES 1 AND 4. VALUE AMPLIFIED BY CONSTANT Z. VALUE IN AMPERES.\n");
fprintf(file4, "Gb 6 3 (2,5) %e \n\n",Kb);
fprintf(file4,"*DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X . DEPENDS OF THE CURRENT THATS FLOWS IN THE VAUX SOURCE. IS AMPLIFIED BY CONSTANTE Y. VALUE IN VOLTS.\n");
fprintf(file4, "Hd 5 8 Vaux %e \n\n",Kd);
fprintf(file4,"*INDEPENDENT VOLTAGE SOURCE TO USE IN THE DEPENDENT VOLTAGE SOURCE. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file4, "Vaux 9 7 0 \n\n");
fprintf(file4,"*CAPACITOR. CONNECTED TO NODES X & X. VALUE = 0 VOLTS\n");
fprintf(file4, "C1 6 8 %e \n\n",C);
fprintf(file4,".IC v(6)=%e v(8)=%e\n\n",V(6), V(8));
fprintf(file4, ".END\n\n");
fclose(file4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf ("\n\n-----------Parte 5-----------\n\n");

printf ("voltages_p5_amplitude_TAB\n");
for i=1:length(I)
  Vm(i)=abs(I(i));
  printf("Vm(%d) = %e \n" ,i ,Vm(i));
endfor
printf ("voltages_p5_amplitude_END\n\n");

printf ("voltages_p5_phase_TAB\n");
for i=1:length(I)
  Vph(i)= arg(I(i));
  printf("Vph(%d) = %e \n" ,i ,Vph(i));
endfor
printf ("voltages_p5_phase_END\n");

%time axis: -5 to 0ms with 1us steps
t_neg=-5e-3:(5e-3)/100:0-(5e-3)/100; %s
t= 0:(20e-3)/1000:20e-3;
t_pos=t;
t_total=[t_neg, t_pos];

Vs_pos=1*sin(w*t);
V6_f=Vm(6)*sin(w*t+Vph(6));

Vs_neg=Vs+0*t_neg;
Vs_pos= e.^(-j*(w*t_pos-pi/2));
Vs_total=[Vs_neg, Vs_pos];

V6_neg=V(6)+0*t_neg;
V6_pos = V6_n + V6_f;
V6_total=[V6_neg, V6_pos];

part5 = figure();
plot(t_total*1e3,Vs_total,'b', "linewidth",4);
hold on
plot(t_total*1000, V6_total, "r", "linewidth",4);

xlabel ("t[ms]");
ylabel (" [V]");

h=legend({"V_{s}","V_{6}"});
legend(h,"location", "northeast");
set (gca,"fontsize",16,"linewidth",2);
set (h,"fontsize",16);

print (part5, "part4.eps", "-depsc");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf ("\n\n-----------Parte 6-----------\n\n");

freq=logspace(-1, 6, 200);

for counter=1:size(freq, 2)

w(counter)=2*pi*freq(counter);

Zc = (1/(2*pi*freq(counter)*C*j));

J =[G1 -G1-G2-G3 G2 0 G3 0 0 0;
    0 Kb+G2 -G2 0 -Kb 0 0 0;
    0 Kb 0 0 -Kb-G5 G5+1/Zc 0 -1/Zc;
    0 0 0 -G6 0 0 G6+G7 -G7;
    0 G3 0 G4 -G3-G5-G4 G5+1/Zc G7 -G7-1/Zc;
    1 0 0 -1 0 0 0 0;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1;
    0 0 0 1 0 0 0 0];
    
L=[0; 0; 0; 0; 0; 1; 0; 0];

%J*K=L
K=inverse(J)*L;

Vc_freq(counter) = K(6)-K(8);
V6_freq(counter) = K(6);
Vs_freq(counter) = K(1)-K(4);

endfor

Vc_eixoY_amp = 20*log10(abs(Vc_freq));
V6_eixoY_amp = 20*log10(abs(V6_freq));
Vs_eixoY_amp = 20*log10(abs(Vs_freq));

Vc_eixoY_ang = angle((Vc_freq))*180/pi;
V6_eixoY_ang= angle((V6_freq))*180/pi;
Vs_eixoY_ang = angle((Vs_freq))*180/pi;


part6_amp = figure();
plot (log10(freq), Vc_eixoY_amp, "color",[0.95,0.60,0], "linewidth",4);
hold on
plot (log10(freq), V6_eixoY_amp, "b", "linewidth",4);
plot (log10(freq), Vs_eixoY_amp, "r", "linewidth",4)

xlabel ("log10(frequency) [Hz]");
ylabel ("Amplitude [V]");

h=legend({"V_{c}","V_{6}","V_{s}"});
legend(h,"location", "northeast");
set (gca,"fontsize",16,"linewidth",2);
set (h,"fontsize",16);

axis([-1 6 -90 20]);
print (part6_amp, "part6_amp.eps", "-depsc");
 
part6_ang = figure();
plot (log10(freq), Vc_eixoY_ang, "color",[0.95,0.60,0], "linewidth",4);
hold on
plot (log10(freq), V6_eixoY_ang, "b","linewidth",4);
plot (log10(freq), Vs_eixoY_ang, "r","linewidth",4);

xlabel ("log10(frequency) [Hz]");
ylabel ("Phase [degrees]");

h=legend({"V_{c}","V_{6}","V_{s}"});
legend(h,"location", "northeast");
set (gca,"fontsize",16,"linewidth",2);
set (h,"fontsize",16);

axis([-1 6 -190 15]);
print (part6_ang, "part6_ang.eps", "-depsc");
