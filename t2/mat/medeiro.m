close all
clear all
	
	%%EXAMPLE SYMBOLIC COMPUTATIONS
	
	pkg load symbolic
	pkg load control
	
	syms R1
	syms R2
	syms R3
	syms R4
	syms R5
	syms R6
	syms R7
	
	syms Kb
	syms Kd
	
	syms Vs
	
	syms Vb
	syms Vd
	syms Ib
	syms Id
	
	syms Ic
	syms Vc
	syms C
	
	syms v1
	syms v2
	syms v3
	syms v4
	syms v5
	syms v6
	syms v7
	syms v8
	
	%%
	%%  PASSO 1
	%%
	
	printf("\n\Knots equações:\n");
	(v3-v2)/R2+(v5-v2)/R3+(v1-v2)/R1 == 0
	Ib+(v2-v3)/R2 == 0
	(v5-v6)/R5-Ib-Ic == 0
	(v4-v7)/R6+(v8-v7)/R7 == 0
	printf("\n\SuperKnots equações:\n");
	(v7-v8)/R7-Ic+(v4-v5)/R4+(v2-v5)/R3+(v6-v5)/R5 == 0
	(v7-v4)/R6+(v5-v4)/R4+(v2-v1)/R1+Ic == 0
	printf("\n\Por observação do circuito\n");
	
	v4 == 0
	
	v1-v4 == Vs
	v5-v8 == Vd
	v2-v5 == Vb
	Id = (v4-v7)/R6
	
	v6-v8 == Vc
	Ic == 0
	
	Vd = Kd*Id
	Ib = Kb*Vb
	printf("\n\Como uma das equações do SuperKnots não é necessária, retiramos a maior:\n\n");
	(v7-v8)/R7+(v4-v5)/R4+(v2-v5)/R3+(v6-v5)/R5+Ic== 0
	
	%%  Atribuição de Valores
	printf("\nAtribuição dos valores para o laboratório:\n\n");
	
	R1 = 1.03258022265e3
	R2 = 2.05854281116e3 
	R3 = 3.05658918951e3 
	R4 = 4.12083818633e3 
	R5 = 3.10223748153e3 
	R6 = 2.09909352125e3 
	R7 = 1.01569886691e3 
	
	Vs = 5.19832384287 
	
	Kb = 7.07448059081e-3 
	Kd = 8.22345657857e3 
	
	C = 1.04739543259e-6
	
	Ic = 0
	
	%     v1    , v2              , v3     , v4       , v5        , v6    , v7          , v8    , Vb    , Vd, Ib, Id  , Vc, Ic
	A = [ 0     , 0               , 0     , 1/R6      , 0         , 0     , -1/R7-1/R6  , 1/R7  , 0     , 0 , 0 , 0   , 0 , 0   ;...
	      0     , 0               , 0     , 0         , 1/R5      , -1/R5 , 0           , 0     , 0     , 0 , -1, 0   , 0 , -1  ;...
	      0     , 1/R2            , -1/R2 , 0         , 0         , 0     , 0           , 0     , 0     , 0 , 1 , 0   , 0 , 0   ;...
	      1/R1  , -1/R2-1/R1-1/R3 , 1/R2  , 0         , 1/R3      , 0     , 0           , 0     , 0     , 0 , 0 , 0   , 0 , 0   ;...
	      -1/R1 , 1/R1            , 0     , -1/R4-1/R6, 1/R4      , 0     , 1/R6        , 0     , 0     , 0 , 0 , 0   , 0 , 0   ;...
	      0     , 0               , 0     , 1         , 0         , 0     , 0           , 0     , 0     , 0 , 0 , 0   , 0 , 0   ;...
	      1     , 0               , 0     , -1        , 0         , 0     , 0           , 0     , 0     , 0 , 0 , 0   , 0 , 0   ;...
	      0     , 0               , 0     , 0         , 1         , 0     , 0           , -1    , 0     , -1, 0 , 0   , 0 , 0   ;...
	      0     , 1               , 0     , 0         , -1        , 0     , 0           , 0     , -1    , 0 , 0 , 0   , 0 , 0   ;...
	      0     , 0               , 0     , 1/R6      , 0         , 0     , -1/R6       , 0     , 0     , 0 , 0 , -1  , 0 , 0   ;...
	      0     , 0               , 0     , 0         , 0         , 0     , 0           , 0     , 0     , 1 , 0 , -Kd , 0 , 0   ;...
	      0     , 0               , 0     , 0         , 0         , 0     , 0           , 0     , -Kb   , 0 , 1 , 0   , 0 , 0   ;...
	      0     , 0               , 0     , 0         , 0         , 1     , 0           , -1    , 0     , 0 , 0 , 0   , -1, 0   ;...
	      0     , 0               , 0     , 0         , 0         , 0     , 0           , 0     , 0     , 0 , 0 , 0   , 0 , 1   ]
        
        
        
       B = [ 0; 0; 0; 0; 0; 0; Vs; 0; 0; 0; 0; 0; 0; 0];

x = A\B
v1=x(1)
v2=x(2)
v3=x(3)
v4=x(4)
v5=x(5)
v6=x(6)
v6_Menor0 = v6;
v7=x(7)
v8=x(8)
v8_Menor0 = v8;
Vb=x(9)
Vd=x(10)
Ib=x(11)
Id=x(12)

Vc=x(13)
Ic=x(14)

I1  =(v1-v2)/R1;
I2  =(v2-v3)/R2;
I3  =(v2-v5)/R3;
I4  =(v4-v5)/R4;
I5  =(v5-v6)/R5;
I6  =Id;
I7  =(v7-v8)/R7;
I_S =-I1;
I_D =I7+Ic;

%	Impressão da Tabela
printf("PASSO1_TAB \n")
printf("$v_{1}$ = %e \n", v1);
printf("$v_{2}$ = %e \n", v2);
printf("$v_{3}$ = %e \n", v3);
printf("$v_{4}$ = %e \n", v4);
printf("$v_{5}$ = %e \n", v5);
printf("$v_{6}$ = %e \n", v6);
printf("$v_{7}$ = %e \n", v7);
printf("$v_{8}$ = %e \n", v8);

printf("@$I_{1}$ = %e \n", I1);
printf("@$I_{2}$ = %e \n", I2);
printf("@$I_{3}$ = %e \n", I3);
printf("@$I_{4}$ = %e \n", I4);
printf("@$I_{5}$ = %e \n", I5);
printf("@$I_{6}$ = %e \n", I6);
printf("@$I_{7}$ = %e \n", I7);
printf("@$I_{S}$ = %e \n", I_S);
printf("@$I_{D}$ = %e \n", I_D);

printf("$V_{b}$ = %e \n", Vb);
printf("$V_{d}$ = %e \n", Vd);
printf("@$I_{b}$ = %e \n", Ib);
printf("@$I_{d}$ = %e \n", Id);

printf("$V_{c}$ = %e \n", Vc);
printf("@$I_{c}$ = %e \n", Ic);

printf("PASSO1_END \n")

%%
%%  Passo 2 - Condição limite (t = 0)
%%

% Dizer que Vs = 0
Vs_dummy = 0;

%     v1    , v2              , v3     , v4       , v5        , v6    , v7          , v8    , Vb    , Vd, Ib, Id  , Ic
A = [ 0     , 0               , 0     , 1/R6      , 0         , 0     , -1/R7-1/R6  , 1/R7  , 0     , 0 , 0 , 0   , 0   ;...
      0     , 0               , 0     , 0         , 1/R5      , -1/R5 , 0           , 0     , 0     , 0 , -1, 0   , -1  ;...
      0     , 1/R2            , -1/R2 , 0         , 0         , 0     , 0           , 0     , 0     , 0 , 1 , 0   , 0   ;...
      1/R1  , -1/R2-1/R1-1/R3 , 1/R2  , 0         , 1/R3      , 0     , 0           , 0     , 0     , 0 , 0 , 0   , 0   ;...
      -1/R1 , 1/R1            , 0     , -1/R4-1/R6, 1/R4      , 0     , 1/R6        , 0     , 0     , 0 , 0 , 0   , 0   ;...
      0     , 0               , 0     , 1         , 0         , 0     , 0           , 0     , 0     , 0 , 0 , 0   , 0   ;...
      1     , 0               , 0     , -1        , 0         , 0     , 0           , 0     , 0     , 0 , 0 , 0   , 0   ;...
      0     , 0               , 0     , 0         , 1         , 0     , 0           , -1    , 0     , -1, 0 , 0   , 0   ;...
      0     , 1               , 0     , 0         , -1        , 0     , 0           , 0     , -1    , 0 , 0 , 0   , 0   ;...
      0     , 0               , 0     , 1/R6      , 0         , 0     , -1/R6       , 0     , 0     , 0 , 0 , -1  , 0   ;...
      0     , 0               , 0     , 0         , 0         , 0     , 0           , 0     , 0     , 1 , 0 , -Kd , 0   ;...
      0     , 0               , 0     , 0         , 0         , 0     , 0           , 0     , -Kb   , 0 , 1 , 0   , 0   ;...
      0     , 0               , 0     , 0         , 0         , 1     , 0           , -1    , 0     , 0 , 0 , 0   , 0   ]
      
B = [ 0     ; 0               ; 0     ; 0         ; 0         ; 0     ; Vs_dummy    ; 0     ; 0     ; 0 ; 0 ; 0   ; Vc  ];

x = A\B

Ic = -x(13)

R_eq = Vc/Ic
tau = R_eq*C

v1=x(1)
v2=x(2)
v3=x(3)
v4=x(4)
v5=x(5)
v6=x(6)
v7=x(7)
v8=x(8)
Vb=x(9)
Vd=x(10)
Ib=x(11)
Id=x(12)

printf("PASSO2_TAB \n")

printf("$V_{c}$ = %e \n", Vc);
printf("@$I_{c}$ = %e \n", Ic);

printf("#$R_{eq}$ = %e \n", R_eq);
printf("#$\\tau$ = %e \n", tau);

printf("PASSO2_END \n")

%%
%%  PASSO 3 - Solução Natural - v6
%%

tempo = 0:(20e-3)/1000:20e-3;

v6_Nt = v6*e.^(-1/tau*tempo);

hf_PASSO3 = figure ();

plot(tempo*1e3,v6_Nt)
hold on

xlabel ("t[ms]");
ylabel ("v6(t) [V]");
print (hf_PASSO3, "PASSO4.eps", "-depsc");

%%
%%  PASSO 4 - Solução Forçada - Todos os Nodos (apenas o coeficiente)
%%

Vs_dummy = 1;

syms Zc
syms w

w = 2*pi*1e3
Zc = 1/(j*w*C)

%     v1    , v2              , v3     , v4       , v5        , v6          , v7          , v8    , Vb    , Vd, Ib, Id
A = [ 0     , 0               , 0     , 1/R6      , 0         , 0           , -1/R7-1/R6  , 1/R7  , 0     , 0 , 0 , 0   ;...
      0     , 0               , 0     , 0         , 1/R5      , -1/R5-1/Zc  , 0           , 1/Zc  , 0     , 0 , -1, 0   ;...
      0     , 1/R2            , -1/R2 , 0         , 0         , 0           , 0           , 0     , 0     , 0 , 1 , 0   ;...
      1/R1  , -1/R2-1/R1-1/R3 , 1/R2  , 0         , 1/R3      , 0           , 0           , 0     , 0     , 0 , 0 , 0   ;...
      -1/R1 , 1/R1            , 0     , -1/R4-1/R6, 1/R4      , 0           , 1/R6        , 0     , 0     , 0 , 0 , 0   ;...
      0     , 0               , 0     , 1         , 0         , 0           , 0           , 0     , 0     , 0 , 0 , 0   ;...
      1     , 0               , 0     , -1        , 0         , 0           , 0           , 0     , 0     , 0 , 0 , 0   ;...
      0     , 0               , 0     , 0         , 1         , 0           , 0           , -1    , 0     , -1, 0 , 0   ;...
      0     , 1               , 0     , 0         , -1        , 0           , 0           , 0     , -1    , 0 , 0 , 0   ;...
      0     , 0               , 0     , 1/R6      , 0         , 0           , -1/R6       , 0     , 0     , 0 , 0 , -1  ;...
      0     , 0               , 0     , 0         , 0         , 0           , 0           , 0     , 0     , 1 , 0 , -Kd ;...
      0     , 0               , 0     , 0         , 0         , 0           , 0           , 0     , -Kb   , 0 , 1 , 0   ]
      
B = [ 0     ; 0               ; 0     ; 0         ; 0         ; 0           ; Vs_dummy    ; 0     ; 0     ; 0 ; 0 ; 0   ];

x = A\B


%%  Ex
Vs_Ft     = e.^(j*(w*tempo-pi/2));
v6_Ft     = Vs_Ft*(x(6));

%%
%%  PASSO 5 - Solução Prévia + Natural + Forçada - VS, v6
%%

tempoNegativo = -5e-3:(5e-3)/100:0-(5e-3)/100;
tempoPositivo = tempo;
tempoTotal    = [tempoNegativo, tempoPositivo];
 
vs_negativo = Vs+0*tempoNegativo;
vs_positivo = e.^(j*(w*tempoPositivo-pi/2));
vs_total = [vs_negativo, vs_positivo];

v6_negativo = v6_Menor0+0*tempoNegativo;
v6_positivo = v6_Nt + v6_Ft;
v6_total = [v6_negativo, v6_positivo];

%v8_negativo = v8_Menor0+0*tempoNegativo;
%v8_positivo = v8_Ft;
%v8_total = [v8_negativo, v8_positivo];

hf_PASSO5 = figure ();

plot(tempoTotal*1e3,vs_total,'r')
hold on
plot(tempoTotal*1e3,v6_total,'b')
%plot(tempoTotal,v8_total,'g')
%plot(tempoTotal,v6_total-v8_total, '-')

xlabel ("t[ms]");
ylabel ("vi(t), vo(t) [V]");
legend;
print (hf_PASSO5, "PASSO5.eps", "-depsc");
hold off

%%
%%  PASSO 6 - Resposta de Frequência - Vs, v6, Vc
%%

w = logspace(-1, 6, 200);

Vs_frec_response = w;
Vc_frec_response = w;
v6_frec_response = w;
Vc_frec_response_teste =w;
Vs_dummy = 1;

for c = 1:size(w,2)
  
  Zc = 1/(j*2*pi*w(c)*C);
  
  %     v1    , v2              , v3     , v4       , v5        , v6          , v7          , v8    , Vb    , Vd, Ib, Id
  A = [ 0     , 0               , 0     , 1/R6      , 0         , 0           , -1/R7-1/R6  , 1/R7  , 0     , 0 , 0 , 0   ;...
        0     , 0               , 0     , 0         , 1/R5      , -1/R5-1/Zc  , 0           , 1/Zc  , 0     , 0 , -1, 0   ;...
        0     , 1/R2            , -1/R2 , 0         , 0         , 0           , 0           , 0     , 0     , 0 , 1 , 0   ;...
        1/R1  , -1/R2-1/R1-1/R3 , 1/R2  , 0         , 1/R3      , 0           , 0           , 0     , 0     , 0 , 0 , 0   ;...
        -1/R1 , 1/R1            , 0     , -1/R4-1/R6, 1/R4      , 0           , 1/R6        , 0     , 0     , 0 , 0 , 0   ;...
        0     , 0               , 0     , 1         , 0         , 0           , 0           , 0     , 0     , 0 , 0 , 0   ;...
        1     , 0               , 0     , -1        , 0         , 0           , 0           , 0     , 0     , 0 , 0 , 0   ;...
        0     , 0               , 0     , 0         , 1         , 0           , 0           , -1    , 0     , -1, 0 , 0   ;...
        0     , 1               , 0     , 0         , -1        , 0           , 0           , 0     , -1    , 0 , 0 , 0   ;...
        0     , 0               , 0     , 1/R6      , 0         , 0           , -1/R6       , 0     , 0     , 0 , 0 , -1  ;...
        0     , 0               , 0     , 0         , 0         , 0           , 0           , 0     , 0     , 1 , 0 , -Kd ;...
        0     , 0               , 0     , 0         , 0         , 0           , 0           , 0     , -Kb   , 0 , 1 , 0   ];
        
  B = [ 0     ; 0               ; 0     ; 0         ; 0         ; 0           ; Vs_dummy    ; 0     ; 0     ; 0 ; 0 ; 0   ];

  x = A\B;
  
  Vs_frec_response(c) = x(1)-x(4);
  v6_frec_response(c) = x(6);
  Vc_frec_response(c) = x(6)-x(8);
  %Vc_frec_response_teste(c) = x(5);
  
endfor  
%
w_absissa = log10(w);

Vs_frec_response_amplitude = 20*log10(abs(Vs_frec_response));
v6_frec_response_amplitude = 20*log10(abs(v6_frec_response));
Vc_frec_response_amplitude = 20*log10(abs(Vc_frec_response));

hf_PASSO6 = figure();

plot(w_absissa,Vs_frec_response_amplitude,'r')
hold on
plot(w_absissa,v6_frec_response_amplitude,'b')
plot(w_absissa,Vc_frec_response_amplitude,'g')
%plot(w_absissa,20*log10(Vc_frec_response_teste),'d')

xlabel ("t[ms]");
ylabel ("vi(t), vo(t) [V]");

print (hf_PASSO6, "PASSO6-AMPLITUDE.eps", "-depsc");

Vs_frec_response_angulo = angle(Vs_frec_response)*180/pi;
v6_frec_response_angulo = angle(v6_frec_response)*180/pi;
Vc_frec_response_angulo = angle(Vc_frec_response)*180/pi;
%Vc_frec_response_teste_ang = angle(Vc_frec_response_teste)*180/pi;
hf_PASSO6 = figure();

plot(w_absissa,Vs_frec_response_angulo,'r')
hold on
plot(w_absissa,v6_frec_response_angulo,'b')
plot(w_absissa,Vc_frec_response_angulo,'g')
%plot(w_absissa,Vc_frec_response_teste_ang,'d')

xlabel ("t[ms]");
ylabel ("vi(t), vo(t) [V]");

print (hf_PASSO6, "PASSO6-ANGULO.eps", "-depsc"); 
        
        
        
        
        
        
        
        
        
        
        
        