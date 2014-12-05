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
% Reactions in series in isothermal PFR and CSTR: comparison
%
%              A -> B -> C
%
%-------------------------------------------------------------------------%
close all;
clear all;

% ------------------------------------------------------------------------%
% Initial concentration
% ------------------------------------------------------------------------%
CA0 = 1;    % [uc]

% ------------------------------------------------------------------------%
% Kinetic parameters
% ------------------------------------------------------------------------%
k1  = 2;    % [1/ut]
k2  = 1;    % [1/ut]

% ------------------------------------------------------------------------%
% Optimal values (analytical solutions)
% ------------------------------------------------------------------------%
CSTR_tOpt = 1/sqrt(k1*k2)
PFR_tOpt = log(k2/k1)/(k2-k1)

% ------------------------------------------------------------------------%
% Profiles (analytical solutions)
% ------------------------------------------------------------------------%
np=1000;
t = linspace(0,100*max(CSTR_tOpt, PFR_tOpt),np);

% CSTR
CSTR_CA = CA0./(1+k1*t);
CSTR_CD = CA0*k1*t./(1+k1*t)./(1+k2*t);
CSTR_CU = CA0 - CSTR_CA - CSTR_CD;

% Plug Flow Reactor
PFR_CA = CA0*exp(-k1*t);
PFR_CD = CA0*k1/(k2-k1)*(exp(-k1*t)-exp(-k2*t));
PFR_CU = CA0 - PFR_CA - PFR_CD;

% Conversion
CSTR_XA = (CA0-CSTR_CA)/CA0;
PFR_XA = (CA0-PFR_CA)/CA0;

% Formation rates
CSTR_RA = -k1*CSTR_CA;
CSTR_RD = k1*CSTR_CA - k2*CSTR_CD;
CSTR_RC = k2*CSTR_CD;

PFR_RA = -k1*PFR_CA;
PFR_RD = k1*PFR_CA - k2*PFR_CD;
PFR_RC = k2*PFR_CD;

% Selectivity
CSTR_S = CSTR_RD./(CSTR_RC);
CSTR_Stilde = CSTR_CD./(CSTR_CU);
PFR_S = PFR_RD./(PFR_RC);
PFR_Stilde = PFR_CD./(PFR_CU);

%Yield
CSTR_Y = CSTR_RD./(-CSTR_RA);
CSTR_Ytilde = CSTR_CD./(CA0-CSTR_CA);
PFR_Y = PFR_RD./(-PFR_RA);
PFR_Ytilde = PFR_CD./(CA0-PFR_CA);

% ------------------------------------------------------------------------%
% Plotting solution
% ------------------------------------------------------------------------%

% Concentrations vs time
subplot(2,2,1);
hold all;
plot(t,CSTR_CD,'LineWidth',2);
plot(t,PFR_CD,'LineWidth',2);
title('Concentration of D');
xlabel('time [ut]');
ylabel('concentrations [uc]');
legend('CSTR', 'PFR');
xlim([0 3*max(CSTR_tOpt, PFR_tOpt)]);

% Concentrations vs Conversion
subplot(2,2,2);
hold all;
plot(CSTR_XA,CSTR_CD,'LineWidth',2);
plot(PFR_XA,PFR_CD,'LineWidth',2);
title('Concentration of D');
xlabel('X_A');
ylabel('concentrations [uc]');
legend('CSTR', 'PFR');
xlim([0 1]);

% Yield vs Conversion
subplot(2,2,3);
hold all;
plot(CSTR_XA,CSTR_Ytilde,'LineWidth',2);
plot(PFR_XA,PFR_Ytilde,'LineWidth',2);
title('Overall Yield');
xlabel('X_A');
ylabel('yield');
legend('CSTR', 'PFR');
xlim([0 1]);

% Selectivity vs Conversion
subplot(2,2,4);
hold all;
plot(CSTR_XA(2:np),CSTR_Stilde(2:np),'LineWidth',2);
plot(PFR_XA(2:np),PFR_Stilde(2:np),'LineWidth',2);
title('Overall selectivity');
xlabel('X_A');
ylabel('selectivity');
legend('CSTR', 'PFR');
xlim([0 1]);
