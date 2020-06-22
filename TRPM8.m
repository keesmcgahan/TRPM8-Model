
clc 
clear 
close all

%What variables are
%X=[m,h,n,V]

%Solving the Odes
t_int=[-30:0.1:40];
sol=ode45(@dhodge,t_int,[.05,.6,.3,-65]');


%%generating information for fixed temperature experiments
x=t_int;
TC=zeros(1,length(x))+5;
TCK=zeros(1,length(x))+5+273;



  y=deval(sol,x,1);
  y2=deval(sol,x,2);
  y3=deval(sol,x,3);
  y4=deval(sol,x,4);
   
  %generating information to plot the currents in fixed temperature plots
 Ik=36*y3.^4.*(y4+70); 
 Ina=120*y.^3.*y2.*(y4-50);
 gm8=3;
 dH=-156000;
 dS=-550;
 zm8=.87;
 Em8=0;
 F=96485;
 R=8.315; 
 am8=1./(1+exp((dH-TCK*dS-(zm8*F*y4)/1000)/(R*TCK)));
 Im8=gm8.*am8.*(y4-Em8);
 

 % plotting the fixed temperature simulation
 z=0*x;
   figure(5)
   yyaxis left
   plot(x,y4,'linewidth',2);
   xlabel('Time (ms)','FontSize',12)
   ylabel('Voltage (mV)','FontSize',12)
    axis([0 30 -80 40])
   yyaxis right
   plot(x,Im8,'linewidth',2);
   hold on
   plot(x,z,'--')
   
   ylabel('I_{m8} (\muA/cm^2)', 'FontSize',12)
   axis([0 30 -100 20])
   

 
 %Temp ramp experiment, comment out lines 16-54, uncomment lines 60-69 
 
% TC=zeros(1,80001)+30;
%  for j=1:60000
% TC(1,j)=.01*abs(.1*(j-1)-3000);  
% end
% figure(7)
% plot(x,y4)
% axis([0 8000 -80 40])
% xlabel('Time (ms)','FontSize',12)
% ylabel('Voltage (mV)','FontSize',12)
%  % ylabel('Temperature (\circC)','FontSize',12)

%the odes

function DXdt=dhodge(t,X)

%parameters
Vr=-65;
Ena=Vr+115;
Ek=Vr-12;
El=Vr+10.613;
Eca=145;
gna=120;
gk=36;
gl=.3;
Cm=1;
F=96485;
R=8.315;
Em8=0;

%TempRamp simulation
% if t<=6000
% TC=.01*abs(t-3000);
% TCK=.01*abs(t-3000)+273;
% else
%     TC=30;
%     TCK=30+273;
% end



%Temperature for fixed temperature experiments
TC=5;
TCK=TC+273;
phi=3^((TC-6.3)/10);




%Trpm8 channel params
dH=-156000;
dS=-550;
zm8=.87;
gm8=3;

%variables
m=X(1);
h=X(2);
n=X(3);
V=X(4);

%rateconstants
am=(.1*(Vr-V+25))./(exp((Vr-V+25)/10)-1);
bm=4*exp((Vr-V)/18);
ah=.07*exp((Vr-V)/20);
bh=1./(exp((Vr-V+30)/10)+1);
an=(.01*(Vr-V+10))./(exp((Vr-V+10)/10)-1);
bn=.125*exp((Vr-V)/80);
am8=1./(1+exp((dH-TCK*dS-(zm8*F*V)/1000)/(R*TCK)));


%Currents
Ina=gna*m.^3.*h.*(V-Ena);
Ik=gk*n.^4.*(V-Ek);
Il=gl*(V-El);
IN=0;
Im8=gm8.*am8.*(V-Em8);
%X=[V,m,h,n]

DXdt=[phi*(am.*(1-m)-bm.*m),...
      phi*(ah.*(1-h)-bh.*h),...
      phi*(an.*(1-n)-bn.*n),...
      (1/Cm)*(IN-Ina-Ik-Il-Im8)]';
  
  
  
end