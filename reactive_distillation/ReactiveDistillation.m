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
%   Copyright(C) 2014 Alberto Cuoci                                       |
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
%
% Semi-batch reactor + Reactive distillation: 
%
%              A + B <-> C + D
%
% For details see: Fogler, Chapter 4 (Web modules: reactive distillation)
%
%-------------------------------------------------------------------------%
close all;
clear all;

% ------------------------------------------------------------------------%
% Global data
% ------------------------------------------------------------------------%
global k;
global KC;
global FBin;
global Qin;

% ------------------------------------------------------------------------%
% Definition of functions corresponding to ODE systems to be solved
% ------------------------------------------------------------------------%
myfuns = ReactiveDistillationFunctions;

%*************************************************************************%
% Part 1: pure semi-batch reactor (no reactive distillation)
%*************************************************************************%

% ------------------------------------------------------------------------%
% Data
% ------------------------------------------------------------------------%
FBin = 3;   % B inlet flow rate [mol/min]
CBin = 5;   % inlet concentration of B [mol/dm3]
NA0 = 300;  % initial amount of A [mol]
NB0 = 0;    % initial amount of B [mol]
NC0 = 0;    % initial amount of C [mol]
ND0 = 0;    % initial amount of D [mol]
V0 = 150;   % initial volume [dm3]
dt = 120;   % final time [min]
T = 325;    % temperature [K]

% ------------------------------------------------------------------------%
% Kinetic and thermodynamic data
% ------------------------------------------------------------------------%
k = 8.88e8*exp(-7032.1/T);      % forward kinetic constant [dm3/mol/min]                  
KC = 5.2*exp((-8000/1.978)*((1./298)-(1./T))); % equilibrium constant[-]  

% ------------------------------------------------------------------------%
% Preliminary calculations
% ------------------------------------------------------------------------%
Qin = FBin/CBin;  % volumetric flow rate 

% ------------------------------------------------------------------------%
% Initial conditions
% ------------------------------------------------------------------------%
y0 = [NA0, NB0, NC0, ND0, V0];

% ------------------------------------------------------------------------%
% ODE solution
% ------------------------------------------------------------------------%
[t,y] = ode15s(@myfuns.semibatch,[0 dt], y0);

% ------------------------------------------------------------------------%
% Plotting solution
% ------------------------------------------------------------------------%
figure;
hold all;
plot(t,y(:,1)./y(:,5), 'LineWidth', 2);
plot(t,y(:,2)./y(:,5), 'LineWidth', 2);
plot(t,y(:,3)./y(:,5), 'LineWidth', 2);
plot(t,y(:,4)./y(:,5), 'LineWidth', 2);
title ('Pure semi-batch (no reactive distillation)');
legend('A','B','C','D');
xlabel('time [min]');
ylabel('concentrations [mol/dm3]');

%*************************************************************************%
% Part 2: semi-batch reactor + reactive distillation (only C component)
%*************************************************************************%

global Ptot;
global FAir;
global PvC;
global CCout;

% ------------------------------------------------------------------------%
% Additional data for component C
% ------------------------------------------------------------------------%
TcC=506.5;  % Critical temperature [K]
PcC=4750;   % Critical pressure [kPa]
rhoC = 974; % Liquid density [g/dm3]
MWC = 74;   % C molecular weight [g/mol]

% ------------------------------------------------------------------------%
% Additional data for reactive distillation
% ------------------------------------------------------------------------%
Ptot=101.3; % total pressure [kPa]
FAir=100;   % Air flow rate [mol/min]

% ------------------------------------------------------------------------%
% Preliminary calculations
% ------------------------------------------------------------------------%
CCout = rhoC/MWC;                   % C conc. [mol/dm3]
Tr=T/TcC;                           % Reduced temperature [-]
PvC=PcC*exp(10.703-(11.0088/Tr)-...
    5.4361*(log(Tr))+0.3058*Tr^6); 	% vapor pressure [kPa]

% ------------------------------------------------------------------------%
% Initial conditions
% ------------------------------------------------------------------------%
y0 = [NA0, NB0, NC0, ND0, V0];

% ------------------------------------------------------------------------%
% ODE solution
% ------------------------------------------------------------------------%
[t,y] = ode15s(@myfuns.reactive_distillation_only_C,[0 dt], y0);

% ------------------------------------------------------------------------%
% Plotting solution
% ------------------------------------------------------------------------%
figure;
hold all;
plot(t,y(:,1)./y(:,5), 'LineWidth', 2);
plot(t,y(:,2)./y(:,5), 'LineWidth', 2);
plot(t,y(:,3)./y(:,5), 'LineWidth', 2);
plot(t,y(:,4)./y(:,5), 'LineWidth', 2);
title ('Semi-batch reactor + reactive distillation (only C component)');
legend('A','B','C','D');
xlabel('time [min]');
ylabel('concentrations [mol/dm3]');


%*************************************************************************%
% Part 3: semi-batch reactor + reactive distillation (every component)
%*************************************************************************%

global PvA;
global PvB;
global PvD;

global CAout;
global CBout;
global CDout;

% ------------------------------------------------------------------------%
% Additional data for component A
% ------------------------------------------------------------------------%
TcA=592.7;   % Critical temperature [K]
PcA=5786;    % Critical pressure [kPa]
rhoA = 974; % Liquid density [g/dm3]
MWA = 60.05;   % C molecular weight [g/mol]

% ------------------------------------------------------------------------%
% Additional data for component B
% ------------------------------------------------------------------------%
TcB=512.6;   % Critical temperature [K]
PcB=8092;    % Critical pressure [kPa]
rhoB = 974; % Liquid density [g/dm3]
MWB = 32.04;   % C molecular weight [g/mol]

% ------------------------------------------------------------------------%
% Additional data for component D
% ------------------------------------------------------------------------%
TcD=647.1;   % Critical temperature [K]
PcD=22060;    % Critical pressure [kPa]
rhoD = 974; % Liquid density [g/dm3]
MWD = 18.01;   % C molecular weight [g/mol]

% ------------------------------------------------------------------------%
% Preliminary calculations
% ------------------------------------------------------------------------%
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

% ------------------------------------------------------------------------%
% Initial conditions
% ------------------------------------------------------------------------%
y0 = [NA0, NB0, NC0, ND0, V0];

% ------------------------------------------------------------------------%
% ODE solution
% ------------------------------------------------------------------------%
[t,y] = ode15s(@myfuns.reactive_distillation_every,[0 32], y0);

% ------------------------------------------------------------------------%
% Plotting solution
% ------------------------------------------------------------------------%
figure;
hold all;
plot(t,y(:,1)./y(:,5), 'LineWidth', 2);
plot(t,y(:,2)./y(:,5), 'LineWidth', 2);
plot(t,y(:,3)./y(:,5), 'LineWidth', 2);
plot(t,y(:,4)./y(:,5), 'LineWidth', 2);
title ('Semi-batch reactor + reactive distillation (every component)');
legend('A','B','C','D');
xlabel('time [min]');
ylabel('concentrations [mol/dm3]');

