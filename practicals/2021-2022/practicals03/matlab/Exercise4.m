%-------------------------------------------------------------------------%
%        ___  ___      _   _       _       ___ _____ ______ _____         |
%        |  \/  |     | | | |     | |     /   /  __ \| ___ \  ___|        |
%        | .  . | __ _| |_| | __ _| |__  / /| | /  \/| |_/ / |__          |
%        | |\/| |/ _` | __| |/ _` | '_ \/ /_| | |    |    /|  __|         |
%        | |  | | (_| | |_| | (_| | |_) \___  | \__/\| |\ \| |___         |
%        \_|  |_/\__,_|\__|_|\__,_|_.__/    |_/\____/\_| \_\____/         |
%                                                                         |
%                                                                         |
%   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
%   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
%   Department of Chemistry, Materials and Chemical Engineering           |
%   Politecnico di Milano                                                 |
%   P.zza Leonardo da Vinci 32, 20133 Milano                              |
%                                                                         |
%-------------------------------------------------------------------------|
%                                                                         |
%   This file is part of Matlab4CRE framework.                            |
%                                                                         |
%	License                                                               |
%                                                                         |
%   Copyright(C) 2021 Alberto Cuoci                                       |
%   Matlab4CRE is free software: you can redistribute it and/or modify    |
%   it under the terms of the GNU General Public License as published by  |
%   the Free Software Foundation, either version 3 of the License, or     |
%   (at your option) any later version.                                   |
%                                                                         |
%   Matlab4CRE is distributed in the hope that it will be useful,         |
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
%   GNU General Public License for more details.                          |
%                                                                         |
%   You should have received a copy of the GNU General Public License     |
%   along with Matlab4CRE. If not, see <http://www.gnu.org/licenses/>.    |
%                                                                         |
%-------------------------------------------------------------------------%
%                                                                         %
% Tubular reactor for pyrolysis of propane                                %
% Non-isothermal, with heat exchange                                      %
% Parallel exothermic reactions, constant T external environment          %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all; clear variables;

global A1;          % [1/s]
global Ea1;         % [K]
global deltaHR1;    % [cal/mol]
global A2;          % [1/s]
global Ea2;         % [K]
global deltaHR2;    % [cal/mol]
global P;           % [Pa]
global S;           % [m2]
global CpMix;       % [cal/mol/K]
global mtot;        % [kg/h]
global Di;          % [m]
global De;          % [m]
global lambdaMix;   % [kcal/m/s/K]
global muMix;       % [kg/m/s]
global Tfgas;       % [K]
global he;          % [kcal/m2/s/K]
global lambdaWall;  % [kcal/m/s/K]


% ------------------------------------------------------------------------%
% User data                                                               %
% ------------------------------------------------------------------------%
A1 = 2.2e12;            % frequency factor [1/s]
Ea1 = 29500;            % activation temperature [K]
A2 = 7.3e11;            % frequency factor [1/s]
Ea2 = 28700;            % activation temperature [K]
deltaHR1 = 18750;       % reaction heat [cal/mol]
deltaHR2 = 30840;       % reaction heat [cal/mol]

Di = 0.11;              % internal diameter [m]
De = 0.13;              % external diameter [m]
L  = 70;                % length [m]

% C3H8 C2H4 CH4 C3H6 H2 H2O
MW = [  44.0962 28.0536 16.0426 ...
        42.0804  2.0158 18.0152 ]; % molecular weights of species
    
minlet = [2000 0 0 0 0 1000]/3600; % inlet mass flow rates [kg/s]
Tin = 600+273.15;                  % inlet temperature [K]

P          = 3*101325;          % pressure [Pa]
Tfgas      = 1100+273.15;       % temperature of external gases in the furnace [K]
he         = 0.080;             % external heat exchange coefficient [kcal/m2/s/K]
lambdaWall = 15/3600;           % thermal conductivity of tube wall [kcal/m/s/K]

CpMix = 0.8;            % spcific heat [kcal/kg/K]
muMix = 3.8e-5;         % viscosity [kg/m/s]
lambdaMix = 3.194e-5;   % thermal conductivity [kcal/m/s/K]

% ------------------------------------------------------------------------%
% Calculations                                                            %
% ------------------------------------------------------------------------%

S = pi/4*Di^2;          % cross section area [m2]
Ctotin = P/8314/Tin;    % inlet concentration [kmol/m3]
mtot = sum(minlet);     % total mass flow rate [kg/s]
Win = minlet/mtot;      % inlet mass fractions [-]
MWin = 1/sum(Win./MW);  % inlet molecular weight [kg/kmol]
rhoin = Ctotin*MWin;    % inlet density [kg/m3]
Qtotin = mtot/rhoin;    % inlet volumetric flow rate [m3/s]
Ftotin = Ctotin*Qtotin; % inlet molar flow rate [kmol/s]
Xin = Win*MWin./MW;     % inlet mole fractions [-]
Fin = Ftotin*Xin;       % inlet molar flow rates of single species [kmol/s]

% ------------------------------------------------------------------------%
% ODE solution                                                            %
% ------------------------------------------------------------------------%
Yin = [Fin Tin]';
options = odeset('RelTol',1e-6, 'AbsTol',1e-8);
[z, Y] = ode45(@PFR, [0 L], Yin, options);

% Post processing
F_C3H8 = Y(:,1)*3600;           % molar flow rate [kmol/h]
F_C2H4 = Y(:,2)*3600;           % molar flow rate [kmol/h]
F_CH4  = Y(:,3)*3600;           % molar flow rate [kmol/h]
F_C3H6 = Y(:,4)*3600;           % molar flow rate [kmol/h]
F_H2   = Y(:,5)*3600;           % molar flow rate [kmol/h]
F_H2O  = Y(:,6)*3600;           % molar flow rate [kmol/h]
T      = Y(:,7);                % temperature [K]

% Reconstruction of external wall temperature
for i=1:length(z)
    Ftot = sum(Y(i,1:6));
    [~, Ue] = HeatExchangeCoefficients(Ftot, T(i));
    Twe = Tfgas-Ue/he*(Tfgas-T);      % [K]
end

% Conversion of propane
X = (Fin(1)-F_C3H8/3600)/Fin(1); 

figure;
plot(z,X,'-o');
xlabel('axial coordinate [m]'); ylabel('conversion'); title('conversion');

figure;
plot(z,T-273.15, '-o', z, Twe-273.15, '-o');
xlabel('axial coordinate [m]'); ylabel('T [C]'); 
title('temperature');

figure;
plot(z,F_C3H8,'-o', z, F_C2H4,'-o', z, F_C3H6,'-o', z, F_H2O,'-o');
xlabel('axial coordinate [m]'); ylabel('molar flow rates [kmol/h]'); 
title('molar flow rates'); legend('C3H8', 'C2H4', 'C3H6', 'H2O', 'Location','northeast');

% ------------------------------------------------------------------------%
% ODE system (defined locally, requires at least MATLAB-2016b)            %
% ------------------------------------------------------------------------%
function dY = PFR(z,Y)

    global A1;          % [1/s]
    global Ea1;         % [cal/mol]
    global deltaHR1;    % [cal/mol]
    global A2;          % [1/s]
    global Ea2;         % [cal/mol]
    global deltaHR2;    % [cal/mol]
    global P;           % [Pa]
    global S;           % [m2]
    global CpMix;       % [kcal/kg/K]
    global mtot;        % [kg/s]
    global Di;          % [m]
    global Tfgas;       % [K]

    % Recovering variables
    F_C3H8 = Y(1);   % molar flow rate of propane [kmol/s]
    F_C2H4 = Y(2);   % molar flow rate of ethylene [kmol/s]
    F_CH4  = Y(3);   % molar flow rate of methane [kmol/s]
    F_C3H6 = Y(4);   % molar flow rate of propylene [kmol/s]
    F_H2   = Y(5);   % molar flow rate of hydrogen [kmol/s]
    F_H2O  = Y(6);   % molar flow rate of water steam [kmol/s]
    T      = Y(7);   % temperature [K]
    
    Ftot = F_C3H8+F_C2H4+F_CH4+F_C3H6+F_H2+F_H2O;   % [kmol/s]
    Ctot = P/8314/T;                                % [kmol/m3]
    C_C3H8 = Ctot*F_C3H8/Ftot;                      % [kmol/m3]
    
    % Reaction rates
    k1 = A1*exp(-Ea1/T);            % [1/s];
    k2 = A2*exp(-Ea2/T);            % [1/s];
    r1 = k1*C_C3H8;                 % [kmol/m3/s]
    r2 = k2*C_C3H8;                 % [kmol/m3/s]
    
    % Mass balance equations
    dF_C3H8 = -S*(r1+r2);           % [kmol/m/s]
    dF_C2H4 =  S*r1;                % [kmol/m/s]
    dF_CH4  =  S*r1;                % [kmol/m/s]
    dF_C3H6 =  S*r2;                % [kmol/m/s]
    dF_H2   =  S*r2;                % [kmol/m/s]
    dF_H2O  =  0;                   % [kmol/m/s]
    
    % Heat exchange coefficients
    [Ui, ~] = HeatExchangeCoefficients(Ftot, T);
    
    % Energy balance equation
    QR = -r1*deltaHR1 - r2*deltaHR2;                % [kcal/m3/s]
    dT = -S*( Ui*(T-Tfgas)*4/Di - QR)/(mtot*CpMix); % [K/m]

    % Recover derivatives in time
    dY = [dF_C3H8 dF_C2H4 dF_CH4 dF_C3H6 dF_H2 dF_H2O dT]';
    
end

function [Ui, Ue] = HeatExchangeCoefficients(Ftot, T)

    global P;           % [Pa]
    global De;          % [m]
    global Di;          % [m]
    global S;           % [m2]
    global lambdaMix;   % [kcal/m/s/K]
    global muMix;       % [kg/m/s]
    global CpMix;       % [kcal/kg/K]
    global he;          % [kcal/m2/s/K]
    global lambdaWall;  % [kcal/m/s/K]
    global mtot;        % [kg/s]
    
    Ctot = P/8314/T;
    Qtot = Ftot/Ctot;                   % [m3/s]
    v = Qtot/S;                         % [m/s]
    rho = mtot/(v*S);                   % [kg/m3]
    Re = rho*v*Di/muMix;                % [-]
    Pr = muMix*CpMix/lambdaMix;         % [-]
    Nu = 0.023*Re^0.8*Pr^0.333;         % [-]
    hi = Nu*lambdaMix/Di;               % [kcal/m2/s/K]
	s = (De-Di)/2;                      % [m]
	
    Ue = 1/( De/Di/hi + 1/he + ...
            (s/lambdaWall)*(De*log(De/Di))/(De-Di) );   % [kcal/m2/s/K]
    Ui = Ue*De/Di;                                      % [kcal/m2/s/K]
    
end