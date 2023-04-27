% Parameter estimation by using fmincon funtion

global BW h ib Vg Vm0 Vi Gpb ipb yobs texp
 

texp = [0 60 120 150  180 210  240];
yobs = [70.631 70.631  67.158 317.178 237.311 203.697 155.082];
% amy_exp = [4.736 4.736 4.736 10.752 13.312 12.8 12.288];
%yobs1 = amy_exp;

% constant parameters as a global variable
BW = 78;     %kg Lean
h = 91.76;  %mg/dl  GLUCOSE (Glu_basal)
ib = 67.158;     %25.49; %pmol/l Insulin  basal)
Vg = 1.88;   % Glucose Volume in dl/Kg
Vi = 0.05;   % Insulin Volume in L/kg
Gpb = h*Vg; % mg/dl * dl/Kg = mg/Kg
ipb = ib*Vi; %pmol/kg % amount of plasma insulin


% parameters 
EGPb=1.92;           %parameter(1)
parameter(1) =EGPb; 
kstomach = 0.0558 ;  %parameter(2)
parameter(2) =kstomach;
kgut = 0.160 ;       %parameter(3)
parameter(3) =kgut;
kabs = 0.057 ;       %parameter(4)
parameter(4) =kabs;
f=0.90;              %parameter(5)
parameter(5) =f;
ke1=0.0005;          %parameter(6)
parameter(6) =ke1;
ke2=339;             %parameter(7)
parameter(7) =ke2;
Vmx=0.047;           %parameter(8)
parameter(8) =Vmx;
km0=225.59;          %parameter(9)
parameter(9) =km0;
K2=0.079;            %parameter(10)
parameter(10) =K2;
K1=0.065;            %parameter(11)
parameter(11) =K1;
Fsnc=1;              %parameter(12)
parameter(12) =Fsnc;
p2U=0.0331;          %parameter(13)
parameter(13) =p2U;
K=2.30;              %parameter(14)
parameter(14) =K;
beta=0.021;          %parameter(15)
parameter(15) =beta;
alpha=0.025;         %parameter(16)
parameter(16) =alpha;
gamma=0.39;          %parameter(17)
parameter(17) =gamma;
m1=0.060;            %parameter(18)
parameter(18) =m1;
m2=0.82;             %parameter(19)
parameter(19) =m2;
m5=0.035;            %parameter(20)
parameter(20) =m5;
HEb=0.262;           %parameter(21)
parameter(21) =HEb;
ki=0.0079;	         %parameter(22)           
parameter(22) =ki;
kp2=0.0021;          %parameter(23)        
parameter(23) =kp2;
kp3=0.009;           %parameter(24) 
parameter(24) =kp3;
kp4=0.0618;          %parameter(25)
parameter(25) =kp4;


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

% glucose ingested in grams (MEAL =447 Kcal, 72% glucose -> 4kcal/g -> 80 g
mealamounts = 80 ;                 % glucose ingested in grams
mealamounts = 1000*mealamounts ;     % convert to milligrams

% % Initial conditions
y0=[mealamounts ,0,0,Gpb,Gtb,0,ipo,ilb,ipb,0,ib,ib];

% fmincon
fun = @glucose_insulin_model_lean;
A = [];
b= [];
Aeq = [];
beq = [];

lb =  [1.728, 0.0502,0.1440,0.0513,0.81, 0.00045,305.1,0.0423,...
      203.031,0.0711,0.0585,0.9,0.0298,2.03,0.019,0.023, 0.36,...
      0.055, 0.74,0.035, 0.23, 0.0072, 0.0019, 0.0080, 0.05418];

ub =  [2.11,0.0614,0.176,0.0627,0.99, 0.00055, 372.9, 0.0517,...
       248.15,0.0869, 0.0715,1.1,0.0364,2.53,0.023,0.028,0.42,...
       0.066, 0.90, 0.037, 0.28, 0.0083, 0.0023, 0.0095, 0.07];


nonlcon = @unitdisk;
x0 = parameter;
% options = optimset('Algorithm','interior-point');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% options = optimoptions('fmincon',...
%     'Algorithm','sqp','Display','iter','ConstraintTolerance',1e-13);
% options = optimoptions(options,'StepTolerance',1e-13);
options = optimset('PlotFcns','optimplotfval','TolX',1e-13);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[param, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);



    
