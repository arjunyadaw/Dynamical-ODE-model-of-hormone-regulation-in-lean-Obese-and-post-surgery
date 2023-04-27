function derivative=dydt_lean(t,y,parameter)

global BW h ib Vg Vm0 Vi Gpb ipb yobs texp

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

   

m4=2/5*m2*HEb; %min^-1
ilb=ipb*(m4+m2)/m1; %pmol/kg
m3=HEb*m1/(1-HEb);
SRb=ipb*m4+ilb*m3; %pmol/kg/min
ipo=SRb/gamma; %pmol/kg
m6=m5*SRb+HEb;

kp1=EGPb+kp2*Gpb+kp3*ib+kp4*ipo;

G_stomach_1 = y(1) ;
G_stomach_2 = y(2) ;
G_intestine = y(3) ;
G_plasma = y(4) ;
G_tissues = y(5) ;
I_free = y(6) ;
I_portal = y(7) ;
I_liver = y(8) ;
I_plasma = y(9) ;
I_interstitial = y(10) ;
I_delay_1 = y(11) ;
I_delay_2 = y(12) ;

% % % % % % % % % % % % % % % % % 
% GLUCOSE
% % % % % % % % % % % % % % % % % 

% Absorption -> Stomach -> Intestine

dG_stomach_1 = -kstomach*G_stomach_1; 
dG_stomach_2=kstomach*G_stomach_1-G_stomach_2*kgut; 
dG_intestine=kgut*G_stomach_2-kabs*G_intestine; 

% PLASMA -> KIDNEY EXCRETION
if (G_plasma > ke2)  % ke2 = threshold for excretion  
  Excretion=ke1*(G_plasma-ke2) ; 
else
  Excretion = 0 ;
end

% INTESTINE -> PLASMA
dG_plasma = max((kp1-kp2*G_plasma-kp3*I_delay_2-kp4*I_portal),0) + ...
            max(kabs*G_intestine*f/(BW),0)+...
            -Fsnc-K1*G_plasma+K2*G_tissues-Excretion ;
% Glucose in tissues
dG_tissues = K1*G_plasma-K2*G_tissues - ...
             (Vm0+Vmx*I_interstitial)/(km0+G_tissues)*G_tissues;

%% I_free = insulin in plasma
if ( (beta*(G_plasma/Vg-h))>-SRb )
  %if glucose levels are high enough to trigger beta cells  
  dI_free = -alpha*I_free + alpha*beta*(G_plasma/Vg-h) ; 
  % insuline production in beta cells proportional to plasma glucose
else
  %otherwise they make basal insulin  
  dI_free=-alpha*I_free-alpha*SRb; 
end

% if plasma glucose levels are above basal levels 
if (dG_plasma>0)&&(G_plasma/Vg>h)  
  dI_portal = I_free-gamma*I_portal + K*dG_plasma/Vg+SRb ;
  % insulin into portal vein depends on glucose rate of change
else
  dI_portal = I_free - gamma*I_portal + SRb ; 
end

% basal hepatic insulin extraction rate
HE = -m5*gamma*I_portal+m6; 
% liver extracts a smaller fraction when secretion is higher
m3=m1*HE/(1-HE); 

% rate constant of insulin degradation, highest when hepatic extraction is high
dI_liver = gamma*I_portal - (m1+m3)*I_liver + m2*I_plasma; 

dI_plasma = m1*I_liver - (m2+m4)*I_plasma; 

dI_interstitial = -p2U*I_interstitial + p2U*(I_plasma/Vi-ib); 
    
dI_delay_1=ki*(I_plasma/Vi-I_delay_1); 
dI_delay_2=ki*(I_delay_1-I_delay_2); 



% Put all derivatives together into a single column vector

derivative = [ ...
 dG_stomach_1; ...
 dG_stomach_2; ...
 dG_intestine; ...
 dG_plasma; ...
 dG_tissues; ...
 dI_free; ...
 dI_portal; ...
 dI_liver; ...
 dI_plasma; ...
 dI_interstitial; ...
 dI_delay_1; ...
 dI_delay_2] ;

end