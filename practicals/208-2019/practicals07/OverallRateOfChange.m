function [r,E,MH,resL,resG,resK] = OverallRateOfChange(pA, CB)

b=2;                        % [-]
GammaA=1e-6/3600;           % [m2/s]
GammaB=1e-6/3600;           % [m2/s]
fL=0.98;                    % [-]
HA=1e5*1e3;                 % [Pa.m3/kmol]
a=20;                       % [-]
KLa=20/3600.;               % [1/s]
KGa=0.01*(1e-3/3600);       % [kmol/m3/s/Pa]
k=100e6/1e-6/3600.;         % [m6/kmol^2/s]

kPrime=k*CB;                    % [m3/kmol/s]
KL = KLa/a;                     % [m/s]
MH=sqrt(kPrime*GammaA*CB)/KL;   % [-]
resG=1/(KGa);                   % [m3.s.Pa/kmol]
resK=HA/(kPrime*CB*fL);         % [m3.s.Pa/kmol]

if (MH<1)               % Regions I and II
    if (MH<0.3)
        E = 1;          % [-]
    else
        E = 1+MH^2/3;   % [-]
    end
    
    resL = HA/KLa/E;                % [m3.s.Pa/kmol]
    r = pA/(resG+resL+resK);        % [kmol/m3/s]
    
else % Regions III
    
    pAI = pA;                       % [Pa]
    for i=1:100
        Ei = 1+GammaB/GammaA*CB*HA/pAI/b;   % [-]
        
        if (MH<5*Ei)
            E = MH;
        else
            E = Ei;
        end
        
        resL = HA/KLa/E;            % [m3.s.Pa/kmol]
        r = pA/(resG+resL+resK);    % [kmol/m3/s]
        
        rGas = KGa*(pA-pAI);        % [kmol/m3/s]
        if ( abs(r-rGas)/r<0.001)
            break;
        end
        
        pAI = pA-r/(KGa);           % [Pa]
    end 
end
