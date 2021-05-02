close all
clear all

pkg load symbolic;
format long;

%chosen variable values
R1 = 1300
R2 = 5300
C = 0.0005
%given variables
freq=50;
Vp = 230;
%primary/secondary circuit
n = 10;
Vs = Vp/n;
%voltage regulator 
num_diodes = 22;
von = 0.6;
I_s = 1e-14;
Vt = 0.025;
mat_sil = 1;

T = 1/(2*freq);
t_off = (1/4)*T;
w = 2*pi*freq;

for i = 1:20
  f = (Vp/n)*C*w*sin(w*t_off) - (1/R1)*(Vp/n)*cos(w*t_off) - I_s*(exp(12/(mat_sil*Vt*num_diodes))-1);
  fl = (Vp/n)*C*(w^2)*cos(w*t_off)+(1/R1)*(Vp/n)*w*sin(w*t_off);
  t_off = t_off - (f/fl);
end

t_on = (3/4)*T;

Req = 1/((1/R1)+(1/R2));

for i = 1:20
  f = (Vp/n)*cos(w*t_on)+(Vp/n)*cos(w*t_off)*exp(-(1/(Req*C))*(t_on-t_off));
  fl = -w*(Vp/n)*sin(w*t_on)-(Vp/n)*cos(w*t_off)*(1/(Req*C))*exp(-(1/(Req*C))*(t_on-t_off));
  t_on = t_on - (f/fl);
end

t = 0:(1e-6):0.2;

l = length(t);

vOenv = ones(1,l);

for i = 1:l
  if t(i)<=t_off
    vOenv(i) = abs((230/n)*cos(w*t(i)));
  else
    if t(i)<=t_on
      vOenv(i) = (230/n)*abs(cos(w*t_off))*exp(-(1/(Req*C))*(t(i)-t_off));
    else
      t_off = t_off + T;
      t_on = t_on + T;
      if t(i)<=t_off
        vOenv(i) = abs((230/n)*cos(w*t(i)));
      else
        if t(i)<=t_on
          vOenv(i) = (230/n)*abs(cos(w*t_off))*exp(-(1/(Req*C))*(t(i)-t_off));
        end
      end
    end
  end    
end


dc_vOenv = mean(vOenv)
ripple_env = max(vOenv) - min(vOenv);
vOenv_centro = (ripple_env/2) + min(vOenv);

vr = zeros(1, length(t));
vs = Vs * cos(w*t);
for i=1:length(t)
	vr(i) = abs(vs(i));
endfor

rd = (mat_sil*Vt)/(I_s*exp((12/num_diodes)/(mat_sil*Vt)));
ac_vOreg = ((num_diodes*rd)/(num_diodes*rd + R2))*(vOenv - dc_vOenv);

if vOenv_centro >= 12
    dc_vOreg = 12;
else
    dc_vOreg = vOenv_centro;
end

vOreg = ac_vOreg + dc_vOreg;


printf ("values_chosen_TAB\n");
printf ("R1 = %e \n", R1);
printf ("R2 = %e \n", R2);
printf ("C = %e \n", C);
printf ("values_chosen_END\n\n");

ripple_env = max(vOenv) - min(vOenv)
average_env = mean(vOenv)

printf ("envelope_TAB\n");
printf ("RippleEnvelope = %e \n", ripple_env);
printf ("AverageEnvelope = %e \n", average_env);
printf ("envelope_END\n\n");

average_reg = mean(vOreg)
ripple_reg = max(vOreg)-min(vOreg)

printf ("regulator_TAB\n");
printf ("RippleRegulator = %e \n", ripple_reg);
printf ("AverageRegulator= %e \n", average_reg);
printf ("regulator_END\n\n");


%plots
	
%output voltages at rectifier, envelope detector and regulator
hfc = figure(1);
title("Regulator and envelope output voltage v_o(t)")
plot(t*1000, vr,"linewidth",4, ";vo_{rectifier}(t);", t*1000,vOenv, "linewidth",4,";vo_{envelope}(t);", t*1000,vOreg,"color","g","linewidth",4,";vo_{regulator}(t);");
xlabel ("t[ms]")
ylabel ("v_O [Volts]")
legend('Location','northeast');
print (hfc, "all_vout.eps", "-depsc");
	
%Deviations (vOenv - 12) 
hfc = figure(2);
title("Deviations from desired DC voltage")
plot (t*1000,vOenv-12, "linewidth",4, ";v4-12 (t);", t*1000,vOreg-12,"linewidth",4, ";v5-12 (t);");
xlabel ("t[ms]")
ylabel ("v_O [Volts]")
legend('Location','northeast');
print (hfc, "deviation.eps", "-depsc");

cost = (R1+R2)/1000 + C*(10^6) + 0.1*(num_diodes + 4)
MERIT= 1/ (cost* (ripple_reg + abs(average_reg-12) + 10e-6))

printf ("merit_TAB\n");
printf ("Cost = %e \n", cost);
printf ("Merit= %e \n", MERIT);
printf ("merit_END\n\n");





