function [error, ypred_ins] = glucose_insulin_model_lean(parameter)

global BW h ib Vg Vm0 Vi Gpb ipb yobs texp

texp = [0 60 120 150  180 210  240];
yobs = [70.631 70.631  67.158 317.178 237.311 203.697 155.082];


% constant parameters as a global variable
BW = 78;     %kg Lean
h = 91.76;  %mg/dl  GLUCOSE (Glu_basal)
ib = 67.158;     %25.49; %pmol/l Insulin  basal)
Vg = 1.88;   % Glucose Volume in dl/Kg
Vi = 0.05;   % Insulin Volume in L/kg
Gpb = h*Vg; % mg/dl * dl/Kg = mg/Kg
ipb = ib*Vi; %pmol/kg % amount of plasma insulin

% parameters 

EGPb = parameter(1);
kstomach = parameter(2);
kgut = parameter(3);
kabs = parameter(4);
f=parameter(5);
ke1=parameter(6);
ke2=parameter(7);
Vmx=parameter(8);
km0=parameter(9);
K2=parameter(10);
K1=parameter(11);
Fsnc=parameter(12);
p2U=parameter(13);
K=parameter(14);
beta=parameter(15);
alpha=parameter(16);
gamma=parameter(17);
m1=parameter(18);
m2=parameter(19);
m5=parameter(20);
HEb=parameter(21);
ki=parameter(22);       
kp2=parameter(23);      
kp3=parameter(24);
kp4=parameter(25);


% TO CALCULATE BASAL LEVEL!!!
if Gpb<=ke2
    % steady state plasma glucose levels are below the
    % threshold for renal excretion
    Gtb=(Fsnc-EGPb+K1*Gpb)/K2; %mg/kg
    Vm0=(EGPb-Fsnc)*(km0+Gtb)/Gtb; %mg/kg/min
    Rdb=EGPb; %mg/kg/min
    PCRb=Rdb/h; %dl/kg/min
else
    % excretion takes place
    Gtb=((Fsnc-EGPb+ke1*(Gpb-ke2))/Vg+K1*Gpb)/K2;%mg/kg
    Vm0=(EGPb-Fsnc-ke1*(Gpb-ke2))*(km0+Gtb)/Gtb; %mg/kg/min
    Rdb=EGPb-ke1*(Gpb-ke2); %mg/kg/min
    PCRb=Rdb/h; %dl/kg/min
end

m4=2/5*m2*HEb; %min^-1
ilb=ipb*(m4+m2)/m1; %pmol/kg
m3=HEb*m1/(1-HEb);
SRb=ipb*m4+ilb*m3; %pmol/kg/min
ipo=SRb/gamma; %pmol/kg
m6=m5*SRb+HEb;

kp1=EGPb+kp2*Gpb+kp3*ib+kp4*ipo;

tsim = 240 ;
mealtimes = 120 ;           % convert to minutes
% glucose ingested in grams (MEAL =447 Kcal, 72% glucose -> 4kcal/g -> 80 g
mealamounts = 80;                 
mealamounts = 1000*mealamounts ;     % convert to milligrams

boluses = [0,mealamounts] ;

intervals = [[0,mealtimes]',[mealtimes,tsim]'] ;
simulationintervals = length(mealtimes) + 1 ;

% % Initial conditions
y0=[mealamounts ,0,0,Gpb,Gtb,0,ipo,ilb,ipb,0,ib,ib];

% These will be vectors that hold all the results
T = 0 ;
Y = y0 ;

% tspan = texp;

for j=1:simulationintervals 
  
  % % With each meal, add bolus defined by 'mealamounts' to stomach
  % % compartment
  y0(1) = y0(1) + boluses(j) ; 
  % Numerical integration of the model, simulation of dydt
[tempT,tempY]=ode45(@(t,y)dydt_lean(t,y,parameter),intervals(j,:),y0); 
  % End of this interval = new initial conditions
  y0 = tempY(end,:) ;
  % Concatenate to keep track of all time points
  T = [T;tempT(2:end)];
  Y = [Y;tempY(2:end,:)] ;
%   Y = vertcat(Y,tempY(2:end,:));

end

I_plasma = pchip(T,Y(:,9)/Vi,texp);

ypred_ins = I_plasma;

error = sum((ypred_ins-yobs).^2);

end
