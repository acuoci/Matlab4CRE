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
%   Copyright(C) 2017 Alberto Cuoci                                       |
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
% Non isothermal, adiabatic plug flow reactor (constant density)          %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%

close all;
clear variables;

% ------------------------------------------------------------------------%
% Data                                                                    %
% ------------------------------------------------------------------------%
Tin = 24+273.15; % [K]
T0  = 20+273.15; % [K]
Tcr = 53+273.15; % [K]

V = 1135;        % [l]

FA0 = 19.50;     % [kmol/h]
FB0 = 364;       % [kmol/h]
FC0 = 0;         % [kmol/h]
FM0 = 32.60;     % [kmol/h]

QA0 = 1320;      % [m3/h]
QB0 = 6600;      % [m3/h]
QC0 = 0;         % [m3/h]
QM0 = 1320;      % [m3/h]

A = 16.96e12;    % [1/h]
E = 72000;       % [J/mol]

CpA = 146;       % [J/mol/K]
CpB = 75;        % [J/mol/K]
CpC = 192;       % [J/mol/K]
CpM = 82;        % [J/mol/K]

HA0 = -148918;   % [J/mol]
HB0 = -275000;   % [J/mol]
HC0 = -505360;   % [J/mol]

% ------------------------------------------------------------------------%
% Calculations                                                            %
% ------------------------------------------------------------------------%

Qtot = QA0+QB0+QC0+QM0; % [m3/s]
CA0 = FA0/Qtot;         % [kmol/m3]
CB0 = FB0/Qtot;         % [kmol/m3]
CC0 = FC0/Qtot;         % [kmol/m3]
CM0 = FM0/Qtot;         % [kmol/m3]

TetaA0 = CA0/CA0;       % [-]
TetaB0 = CB0/CA0;       % [-]
TetaC0 = CC0/CA0;       % [-]
TetaM0 = CM0/CA0;       % [-]

CpIn = TetaA0*CpA+TetaB0*CpB+TetaC0*CpC+TetaM0*CpM;  % [J/mol/K]
deltaCp = CpC-CpA-CpB;                               % [J/mol/K]
deltaHR0 = HC0-HA0-HB0;                              % [J/mol]

Tau = V/Qtot;   % [h]


% ------------------------------------------------------------------------%
% Graphical solution                                              %
% ------------------------------------------------------------------------%

T = 40:100;
T = T +273.15;
Xmb = A*exp(-E/8.314./T)*Tau./(1+A*exp(-E/8.314./T)*Tau);
Xeb = -CpIn*(T-Tin)./(deltaHR0+deltaCp*(T-T0));

figure;
plot((T-273.15),Xmb,'-o', (T-273.15), Xeb,'-o');
legend('Xmb', 'Xeb');
xlabel('T [C]');
ylabel('X');
title('Main results');
