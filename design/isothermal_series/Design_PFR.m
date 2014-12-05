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
% Reactions in series in isothermal PFR
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

% Optimal residence time
tOpt = log(k2/k1)/(k2-k1);

% Concentrations
CAOpt = CA0*exp(-k1*tOpt);
CDOpt = CA0*k1/(k2-k1)*(exp(-k1*tOpt)-exp(-k2*tOpt));
CUOpt = CA0 - CAOpt - CDOpt;

% Conversion
XOpt = (CA0-CAOpt)/CA0;

% Formation rates
RAOpt = -k1*CAOpt;
RDOpt = k1*CAOpt - k2*CDOpt;
RCOpt = k2*CDOpt;

% Selectivity
SOpt = RDOpt/RCOpt
StildeOpt = CDOpt/CUOpt;

% Fractional Yield
YOpt = RDOpt/(-RAOpt);
YtildeOpt = CDOpt/(CA0-CAOpt);

% ------------------------------------------------------------------------%
% Profiles (analytical solutions)
% ------------------------------------------------------------------------%

np=200;
t = linspace(0,5*tOpt,np);

CA = CA0*exp(-k1*t);
CD = CA0*k1/(k2-k1)*(exp(-k1*t)-exp(-k2*t));
CU = CA0 - CA - CD;

XA = (CA0-CA)/CA0;

RA = -k1*CA;
RD = k1*CA - k2*CD;
RC = k2*CD;

S = RD./(RC+1e-6);
Stilde = CD./(CU+1e-6);

Y = RD./(-RA);
Ytilde = CD./(CA0-CA);


% ------------------------------------------------------------------------%
% Plotting solution
% ------------------------------------------------------------------------%

% Concentrations vs time
subplot(2,2,1);
hold all;
plot(t,CA,'LineWidth',2);
plot(t,CD,'LineWidth',2); 
plot(t,CU,'LineWidth',2);
title('Concentrations');
xlabel('time [ut]');
ylabel('concentrations [uc]');
legend('A', 'D', 'U');
xlim([0 3*tOpt]);

% Concentrations vs Conversion
subplot(2,2,2);
hold all;
plot(XA,CA,'LineWidth',2);
plot(XA,CD,'LineWidth',2);
plot(XA,CU,'LineWidth',2);
title('Concentrations');
xlabel('X_A');
ylabel('concentrations [uc]');
legend('A', 'D', 'U');

% Fractional yield vs Conversion
subplot(2,2,3);
hold all;
plot(XA,Ytilde,'LineWidth',2);
title('Yield');
xlabel('X_A');
ylabel('yield');

% Selectivity vs Conversion
subplot(2,2,4);
hold all;
plot(XA(2:np),S(2:np), 'LineWidth',2);
plot(XA(2:np), Stilde(2:np), 'LineWidth',2);
title('Selectivity');
xlabel('X_A');
ylabel('selectivity');
legend('instantaneous', 'overall');