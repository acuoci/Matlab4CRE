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
%   License                                                               |
%                                                                         |
%   Copyright(C) 2019 Alberto Cuoci                                       |
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
% Semi-batch reactor + Reactive distillation:                             %
%                                                                         %
%              A + B <-> C + D                                            %
%                                                                         %
% For details see: Fogler, Chapter 4 (Web modules: reactive distillation) %
%                                                                         %
%-------------------------------------------------------------------------%
close all; clear variables;


%% Case 1: pure semi-batch reactor (no reactive distillation)

% Operating conditions
FBin = 3;   % B inlet flow rate [mol/min]
CBin = 5;   % inlet concentration of B [mol/dm3]
NA0 = 300;  % initial amount of A [mol]
NB0 = 0;    % initial amount of B [mol]
NC0 = 0;    % initial amount of C [mol]
ND0 = 0;    % initial amount of D [mol]
V0 = 150;   % initial volume [dm3]
dt = 120;   % final time [min]
T = 325;    % temperature [K]

% Kinetic and thermodynamic data
k = 8.88e8*exp(-7032.1/T);      % forward kinetic constant [dm3/mol/min]                  
KC = 5.2*exp((-8000/1.978)*((1./298)-(1./T))); % equilibrium constant[-]  

% Preliminary calculations
Qin = FBin/CBin;  % volumetric flow rate 

% Initial conditions
y0 = [NA0, NB0, NC0, ND0, V0];

% ODE solution
[t,y] = ode15s(@semibatch,[0 dt], y0, [], k,KC,FBin,Qin);

% Plotting solution
figure;
hold all;
plot(t,y(:,1)./y(:,5), 'LineWidth', 2);
plot(t,y(:,2)./y(:,5), 'LineWidth', 2);
plot(t,y(:,3)./y(:,5), 'LineWidth', 2);
plot(t,y(:,4)./y(:,5), 'LineWidth', 2);
title ('Pure semi-batch (no reactive distillation)');   legend('A','B','C','D');
xlabel('time [min]');   ylabel('concentrations [mol/dm^3]');


%% Case 2: semi-batch reactor + reactive distillation (only C component)

% Additional data for component C
TcC=506.5;  % Critical temperature [K]
PcC=4750;   % Critical pressure [kPa]
rhoC = 974; % Liquid density [g/dm3]
MWC = 74;   % C molecular weight [g/mol]

% Additional data for reactive distillation
Ptot=101.3; % total pressure [kPa]
FAir=100;   % Air flow rate [mol/min]

% Preliminary calculations
CCout = rhoC/MWC;                   % C conc. [mol/dm3]
Tr=T/TcC;                           % Reduced temperature [-]
PvC=PcC*exp(10.703-(11.0088/Tr)-...
    5.4361*(log(Tr))+0.3058*Tr^6); 	% vapor pressure [kPa]

% Initial conditions
y0 = [NA0, NB0, NC0, ND0, V0];

% ODE solution
[t,y] = ode15s(@reactive_distillation_only_C,[0 dt], y0, [], ...
                k,KC,FBin,Qin,Ptot,FAir,PvC,CCout);

% Plotting solution
figure;
hold all;
plot(t,y(:,1)./y(:,5), 'LineWidth', 2);
plot(t,y(:,2)./y(:,5), 'LineWidth', 2);
plot(t,y(:,3)./y(:,5), 'LineWidth', 2);
plot(t,y(:,4)./y(:,5), 'LineWidth', 2);
title ('Semi-batch reactor + reactive distillation (only C component)');
legend('A','B','C','D');
xlabel('time [min]');   ylabel('concentrations [mol/dm^3]');


%% Case 3: semi-batch reactor + reactive distillation (every component)

% Additional data for component A
TcA=592.7;      % Critical temperature [K]
PcA=5786;       % Critical pressure [kPa]
rhoA = 974;     % Liquid density [g/dm3]
MWA = 60.05;    % C molecular weight [g/mol]

% Additional data for component B
TcB=512.6;      % Critical temperature [K]
PcB=8092;       % Critical pressure [kPa]
rhoB = 974;     % Liquid density [g/dm3]
MWB = 32.04;    % C molecular weight [g/mol]

% Additional data for component D
TcD=647.1;      % Critical temperature [K]
PcD=22060;      % Critical pressure [kPa]
rhoD = 974;     % Liquid density [g/dm3]
MWD = 18.01;    % C molecular weight [g/mol]

% Preliminary calculations
CAout = rhoA/MWA;                   % A conc. [mol/dm3]
CBout = rhoB/MWB;                   % B conc. [mol/dm3]
CDout = rhoD/MWD;                   % D conc. [mol/dm3]

Tr=T/TcA;                            % Reduced temperature (A) [-]
PvA=PcA*exp(12.446-(12.8016/Tr)-...
    7.1135*(log(Tr))+0.3556*Tr^6); 	 % vapor pressure (A) [kPa]

Tr=T/TcB;                            % Reduced temperature (B) [-]
PvB=PcB*exp(14.413-(14.8248/Tr)-...
    8.8032*(log(Tr))+0.4118*Tr^6); 	 % vapor pressure (B) [kPa]

Tr=T/TcD;                            % Reduced temperature (D) [-]
PvD=PcD*exp(11.060-(11.3760/Tr)-...
    5.9233*(log(Tr))+0.316*Tr^6); 	 % vapor pressure (D) [kPa]

% Initial conditions
y0 = [NA0, NB0, NC0, ND0, V0];

% ODE solution
[t,y] = ode15s(@reactive_distillation_every,[0 dt/4], y0, [], ...
                k,KC,FBin,Qin,Ptot,FAir,[PvA PvB PvC PvD],[CAout CBout CCout CDout]);

% Plotting solution
figure;
hold all;
plot(t,y(:,1)./y(:,5), 'LineWidth', 2);
plot(t,y(:,2)./y(:,5), 'LineWidth', 2);
plot(t,y(:,3)./y(:,5), 'LineWidth', 2);
plot(t,y(:,4)./y(:,5), 'LineWidth', 2);
title ('Semi-batch reactor + reactive distillation (every component)');
legend('A','B','C','D');
xlabel('time [min]');   ylabel('concentrations [mol/dm^3]');



%% ODE system 1: pure semi-batch reactor (no reactive distillation)
function dy = semibatch(t,y, k,KC,FBin,Qin)
    
    V = y(5);                   % volume of liquid phase [dm3]
    CA = y(1)/V;                % concentration of species A [mol/dm3]
    CB = y(2)/V;                % concentration of species B [mol/dm3]
    CC = y(3)/V;                % concentration of species C [mol/dm3]
    CD = y(4)/V;                % concentration of species D [mol/dm3]
    
    R  = k*(CA*CB-CC*CD/KC);    % reaction rate [mol/dm3/min]
    
    dy = zeros(5,1);
    dy(1) = -R*V;               % dNA_over_dt [mol/min]
    dy(2) = -R*V + FBin;        % dNB_over_dt [mol/min]
    dy(3) = R*V;                % dNC_over_dt [mol/min]
    dy(4) = R*V;                % dND_over_dt [mol/min]
    dy(5) = Qin;                % dV_over_dt  [dm3/min]
    
end


%% ODE system 2: semi-batch reactor + reactive distillation (only C component)
function dy = reactive_distillation_only_C(t,y, k,KC,FBin,Qin,Ptot,FAir,PvC,CCout)

    V = y(5);                   % volume of liquid phase [dm3]
    CA = y(1)/V;                % concentration of species A [mol/dm3]
    CB = y(2)/V;                % concentration of species B [mol/dm3]
    CC = y(3)/V;                % concentration of species C [mol/dm3]
    CD = y(4)/V;                % concentration of species D [mol/dm3]
    
    R  = k*(CA*CB-CC*CD/KC);    % reaction rate [mol/dm3/min]
    
    XC = y(3)/(y(1)+y(2)+y(3)+y(4));  % mole fraction of C [-]
    FCout=XC*(PvC/Ptot)*FAir;         % evap. molar flow rate C [mol/min]
    
    dy = zeros(5,1);
    dy(1) = -R*V;                       % dNA_over_dt [mol/min]
    dy(2) = -R*V + FBin;                % dNB_over_dt [mol/min]
    dy(3) =  R*V - FCout;               % dNC_over_dt [mol/min]
    dy(4) =  R*V;                       % dND_over_dt [mol/min]
    dy(5) =  Qin - FCout/CCout;         % dV_over_dt  [dm3/min]
    
end


%% ODE system 3: semi-batch reactor + reactive distillation (every component)
function dy = reactive_distillation_every(t,y, k,KC,FBin,Qin,Ptot,FAir,Pv,Cout)
    
    PvA=Pv(1);      PvB=Pv(2);      PvC=Pv(3);      PvD=Pv(4);
    CAout=Cout(1);  CBout=Cout(2);  CCout=Cout(3);  CDout=Cout(4);

    V = y(5);                   % volume of liquid phase [dm3]
    CA = y(1)/V;                % concentration of species A [mol/dm3]
    CB = y(2)/V;                % concentration of species B [mol/dm3]
    CC = y(3)/V;                % concentration of species C [mol/dm3]
    CD = y(4)/V;                % concentration of species D [mol/dm3]
    
    R  = k*(CA*CB-CC*CD/KC);    % reaction rate [mol/dm3/min]
    
    Ntot = y(1)+y(2)+y(3)+y(4); % total number of moles [mol]
    XA = y(1)/Ntot;             % mole fraction of A [-]
    XB = y(2)/Ntot;             % mole fraction of B [-]
    XC = y(3)/Ntot;             % mole fraction of C [-]
    XD = y(4)/Ntot;             % mole fraction of D [-]
    
    FAout=XA*(PvA/Ptot)*FAir;   % evap. molar flow rate A [mol/min]
    FBout=XB*(PvB/Ptot)*FAir;   % evap. molar flow rate B [mol/min]
    FCout=XC*(PvC/Ptot)*FAir;   % evap. molar flow rate C [mol/min]
    FDout=XD*(PvD/Ptot)*FAir;   % evap. molar flow rate D [mol/min]
    
    dy = zeros(5,1);
    dy(1) = -R*V - FAout;               % dNA_over_dt [mol/min]
    dy(2) = -R*V - FBout + FBin;        % dNB_over_dt [mol/min]
    dy(3) =  R*V - FCout;               % dNC_over_dt [mol/min]
    dy(4) =  R*V - FDout;               % dND_over_dt [mol/min]
    dy(5) =  Qin - FAout/CAout ...
                 - FBout/CBout ...
                 - FCout/CCout ...
                 - FDout/CDout;         % dV_over_dt [dm3/min]
    
end
