%% Function for the calculation or air-sea CO2 flux
function [K,F_CO2, dpCO2]=FCO2(pCO2_aqua, pCO2_atm,T,S,u,Ref)
  
   %INPUTS
   %pCO2_aqua= seawater pCO2 (uatm)
   %pCO2_atm=  atmospheric pCO2 (uatm)
   %T=  Temperature (Celsius)
   %S=  Salinity 
   %u = Wind speed (m/s)
    
   %OUTPUTS
   % F_CO2 =K*a(dpCO2) Air-sea CO2 flux 
   %   
   %K=is the transfer velocity according to the 
   % chosen Ref (See Switch below).
   %dpCO2 is the difference of air and seawater pCO2 
   % a = CO2 solibility constant according to Weiss (1974)
   %%%%% CO2 Transfer velocity calculation %%%%%%%%%       
   
   Sc=Schmidt(T);
   switch Ref
    case 1 % Liss and Merlivat (1986)
        if u<=3.6 
        K = 0.170 .* (Sc./600).^(-2/3) ;    
        elseif u>3.6 & u<=13 
        K = ( 2.85*u-9.65 ) .* ( (Sc./600).^-0.5 );
        else
        K = ( 5.9*u-49.3 ) .* ( (Sc./600).^-0.5 );    
        end
    case 2 % Wanninkhof (1992)
        if u<=6;    
        K=0.31*(u.^2).*((Sc./660).^-0.5); %for slower steadier wind        
        else
        K=0.39*(u.^2).*((Sc./660).^-0.5);
        end
    case 3 % Wanninkhof and McGillis (1999)  short term period ( less than a month)     
        K = 0.0283*u.^3.*( (Sc./660).^-0.5 );
    case 4 % Nightingale et al. (2000)
        K = ( 0.222*u.^2+0.333*u ).*( (Sc./600).^-0.5 );
    case 5 % Jean-Baptiste et Poisson (2002)
        K = 1.45 * u.^1.5 .*( (Sc./310).^-0.5 );
    case 6 % Ho et al. (2006)
        K = 0.266 * u.^2 .*( (Sc./600).^-0.5 );
    case 7 % Sweeney et al. (2007)
        K = ( 0.27*u.^2 ).*( (Sc./660).^-0.5 );
    case 8 % Wanninkhof et al. (2009)
        K = 3+ 0.1*u + 0.064*u.^2 + 0.011*u.^3 .*( (Sc./660).^-0.5 );
    case 9 % Ho et al. (2011) a
        K = 0.286*u.^2 .*( (Sc./600).^-0.5 );
    case 10 % Ho et al. (2011) b
        K = 0.0298*u.^3 .*( (Sc./600).^-0.5 );    
    case 11 % Wanninkhof (2014) 
        K = 0.251*u.^2 .*( (Sc./660).^-0.5 );         
    otherwise
        disp('Wrong Input Values of Reference')
   end
   
   
   dpCO2=pCO2_aqua-pCO2_atm; 
   a=Ko_weiss(T,S);    %Solibility in mmol L^-1 atm^-1 or mmol m^-3 uatm^-1    
   F_CO2 =0.24*K.*a.*dpCO2; %CO2 flux (mmol m^-2 d^-1)
end

%%%%%Subrutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ****Solibuility constant (Weiss, 1974) ***********************************

function [Ko]=Ko_weiss(T,S)
A=[-60.2409, 93.4517, 23.3585];  %mol/Kg.atm
B=[0.023517, -0.023656, 0.0047036]; %mol/Kg.atm
T=T+273.15; %Conversio from Celsius degrees to Kelvins
Ln_Ko=A(1)+(A(2).*(100./T))+(A(3).*log(T./100))+S.*(B(1)+(B(2).*(T./100))+(B(3).*(T./100).^2));
Ko=exp(Ln_Ko);
end
%% ******** Schmidt Number*********
%For water of salinity=35 and temperature range 0-30°C    %%%%%%%%%%%%
function [Sc]=Schmidt(T)
    A = 2073.1;     B = 125.62;     C = 3.6276;     D = 0.043219;
    Sc= A - (B.*T)+(C.*T.^2)-(D.*T.^3);    
end
     
%% From and modified : 
        %   Cecilia Chapa Balcorta 
        %   Ensenada, Baja California
        %   UNIVERSIDAD AUTÓNOMA DE BAJA CALIFORNIA/ UNIVERSIDAD DEL MAR  

%       last version of 11/10/2016 
%                      By BENALLAL MOHAMED ANIS 
%                      IMAGES-ESPACE_DEV, 
%                      UNIVERSITY OF PERPIGNAN VIA DOMITIA
        
