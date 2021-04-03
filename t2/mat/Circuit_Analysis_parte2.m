
pkg load symbolic

% Insertion of the resitors' values generated by datagen.py
R1 = 1.04111259479e+03;
R2 = 2.09945227782e+03;
R3 = 3.13109125645e+03;
R4 = 4.11947040212e+03;
R5 = 3.1155879392e+03;
R6 = 2.04799381798e+03; 
R7 = 1.02754401839e+03;


%dados com numero 95826
%R1 = 1.02591946196e+03; 
%R2 = 2.05138905661e+03;
%R3 = 3.09820234189e+03;
%R4 = 4.07126998238e+03; 
%R5 = 3.04248580013e+03; 
%R6 = 2.00249636786e+03; 
%R7 = 1.0123933318e+03; 
%Vs = 5.22225162486; 
%C = 1.02739994257e-06; 
%Kb = 7.04626963966e-03; 
%Kd = 8.17197962581e+03;

Vs = 5.06871572779;
C = 1.04127523824e-06;
Kb = 7.28747116393e-03;
Kd = 8.11568444746e+03;

% Calculation of the inverse of the resistors. Useful to the Node Method
G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

% Insertion of the remaining values generated by datagen.py


printf ("\n\n-----------Parte 1-----------\n\n");

M =[1 0 0 0 0 0 0 0;
    -G1 G1+G2+G3 -G2 0 -G3 0 0 0;
    0 -Kb-G2 G2 0 Kb 0 0 0;
    0 0 0 1 0 0 0 0;
    0 -G3 0 -G4 G3+G4+G5 -G5 -G7 G7;
    0 Kb 0 0 -Kb-G5 G5 0 0;
    0 0 0 -G6 0 0 G6+G7 -G7;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1];

O=[Vs; 0 ; 0; 0 ; 0 ; 0 ; 0; 0];

% D*F=E <=>
N=inv(M)*O;


% Creating a table in the mat folder with the voltages results of the Node Method
printf ("nos_TAB\n");
printf ("V1 = %e \n", N(1));
printf ("V2 = %e \n", N(2));
printf ("V3 = %e \n", N(3));
printf ("V4 = %e \n", N(4));
printf ("V5 = %e \n", N(5));
printf ("V6 = %e \n", N(6));
printf ("V7 = %e \n", N(7));
printf ("V8 = %e \n", N(8));
printf ("nos_END\n");


printf ("\n\n-----------Parte 2-----------\n\n");

%The capacitor was replaced with a voltage Vx=V6-V8,
%where V6 and V8 are the voltages obtained when t<0
Vx=N(6)-N(8);

% Insertion of the matrices that will allow us to calculate the voltages of 
% each node


A =[1 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0;
    -G1 G1+G2+G3 -G2 0 -G3 0 0 0 0;
    0 -Kb-G2 G2 0 Kb 0 0 0 0;   
    0 -G1 0 0 -G4 0 -G6 0 0;
    0 0 0 0 0 0 G6+G7 -G7 0;
    0 Kb 0 0 -Kb-G5 G5 0 0 1;
    0 0 0 0 0 1 0 -1 0;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1 0];
    
B = [0; 0; 0; 0; 0; 0; 0; Vx; 0];


% A*B=R <=>
R=inverse(A)*B;
%R=A\B;

%computation of the equivalent resistor (Req=Vx/(-Ix))
Req=(Vx/(-R(9)));

%time constant (tau=Req*C)
tau=(Req*C);

% Creating a table in the mat folder with the voltages results
printf ("voltages_p2_TAB\n");
printf ("V1 = %e \n", R(1));
printf ("V2 = %e \n", R(2));
printf ("V3 = %e \n", R(3));
printf ("V4 = %e \n", R(4));
printf ("V5 = %e \n", R(5));
printf ("V6 = %e \n", R(6));
printf ("V7 = %e \n", R(7));
printf ("V8 = %e \n", R(8));
printf ( "voltages_p2_END\n");

printf ("\nnorten_current_p2_TAB\n");
printf ("Ix = %e \n", R(9));
printf ("norton_current_p2_END\n");

printf ("\nReq_tau_p2_TAB\n");
printf ("Req = %e \n", Req);
printf ("tau = %e \n", tau);
printf ("Req_tau_p2_END\n");


%-------------------------------------------------------------------------------------------------------------
printf ("\n\n-----------Parte 3-----------\n\n");


%time axis: 0 to 20ms with 1us steps
t = 0:(20e-3)/1000:20e-3; %s


%Vx is the initial voltage
V6_n=(R(6))*e.^(-t/tau);


%Plot natural solution 
%natural = figure ();
%plot (t*1000, V6_n, "g");
%xlabel ("t[ms]");
%ylabel ("Vnatural(t) [V]");
%print (natural, "natural.eps", "-depsc");

%-------------------------------------------------------------------------------------------------------------
printf ("\n\n-----------Parte 4-----------\n\n");

%frequency [Hz]
f=1000;
w=2*pi*f;
Zc=1/(j*w*C);


D =[G1 -G1-G2-G3 G2 0 G3 0 0 0;
    0 Kb+G2 -G2 0 -Kb 0 0 0;
    0 Kb 0 0 -Kb-G5 G5+1/Zc 0 -1/Zc;
    0 0 0 -G6 0 0 G6+G7 -G7;
    0 G3 0 G4 -G3-G5-G4 G5+1/Zc G7 -G7-1/Zc;
    1 0 0 -1 0 0 0 0;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1;
    0 0 0 1 0 0 0 0];
    
E=[0; 0; 0; 0; 0; 1; 0; 0];


% D*F=E <=>
F=D\E;


% Creating a table in the mat folder with the voltages results
printf ("voltages_p4_TAB\n");
printf ("V1 = %e \n", abs(F(1)));
printf ("V2 = %e \n", abs(F(2)));
printf ("V3 = %e \n", abs(F(3)));
printf ("V4 = %e \n", abs(F(4)));
printf ("V5 = %e \n", abs(F(5)));
printf ("V6 = %e \n", abs(F(6)));
printf ("V7 = %e \n", abs(F(7)));
printf ("V8 = %e \n", abs(F(8)));
printf ( "voltages_p4_END\n");


%-------------------------------------------------------------------------------------------------------------
printf ("\n\n-----------Parte 5-----------\n\n");

printf ("voltages_p5_amplitude_TAB\n");
for i=1:length(F)
  Vm(i)=abs(F(i));
  printf("Vm(%d) = %e \n" ,i ,Vm(i));
endfor
printf ("voltages_p5_amplitude_END\n");


printf ("voltages_p5_phase_TAB\n");
for i=1:length(F)
  Vph(i)= arg(F(i));
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

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!nao e VX e v6 do ponto 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
V6_neg=N(6)+0*t_neg;
V6_pos = V6_n + V6_f;
V6_total=[V6_neg, V6_pos];


part4 = figure();
plot(t_total*1e3,Vs_total,'r');
hold on
plot(t_total*1000, V6_total, "g");

xlabel ("t[ms]");
ylabel (" [V]");
print (part4, "part4.eps", "-depsc");



%-------------------------------------------------------------------------------------------------------------
printf ("\n\n-----------Parte 6-----------\n\n");

freq=logspace(-1, 6, 200);

for counter=1:size(freq, 2)

w(counter)=2*pi*freq(counter);

Zc = (1/(2*pi*freq(counter)*C*j));

G =[G1 -G1-G2-G3 G2 0 G3 0 0 0;
    0 Kb+G2 -G2 0 -Kb 0 0 0;
    0 Kb 0 0 -Kb-G5 G5+1/Zc 0 -1/Zc;
    0 0 0 -G6 0 0 G6+G7 -G7;
    0 G3 0 G4 -G3-G5-G4 G5+1/Zc G7 -G7-1/Zc;
    1 0 0 -1 0 0 0 0;
    0 0 0 Kd*G6 -1 0 -Kd*G6 1;
    0 0 0 1 0 0 0 0];
    
H=[0; 0; 0; 0; 0; 1; 0; 0];

K=inverse(G)*H;
%K=G\H;


Vc_freq(counter) = K(6)-K(8);
V6_freq(counter) = K(6);
Vs_freq(counter) = K(1)-K(4);


w_eixoX=log10(w);


endfor

Vc_eixoY_amp = 20*log10(abs(Vc_freq));
V6_eixoY_amp = 20*log10(abs(V6_freq));
Vs_eixoY_amp = 20*log10(abs(Vs_freq));

Vc_eixoY_ang = angle((Vc_freq))*180/pi;
V6_eixoY_ang= angle((V6_freq))*180/pi;
Vs_eixoY_ang = angle((Vs_freq))*180/pi;


part6_amp = figure();
plot (w_eixoX, Vc_eixoY_amp, "r");
hold on
plot (w_eixoX, V6_eixoY_amp, "g");
plot (w_eixoX, Vs_eixoY_amp, "b")

xlabel ("Vi(t)/Vo(t) [rad/s]");
ylabel ("[dB]");
axis([0 6.5 -90 20]);
print (part6_amp, "part6_amp.eps", "-depsc");
 
 
part6_ang = figure();
plot (w_eixoX, Vc_eixoY_ang, "r");
hold on
plot (w_eixoX, V6_eixoY_ang, "g");
plot (w_eixoX, Vs_eixoY_ang, "b");


xlabel ("Vi(t)/Vo(t) [rad/s]");
ylabel ("[dB]");
axis([0 6.5 -190 15]);
print (part6_ang, "part6_ang.eps", "-depsc");
