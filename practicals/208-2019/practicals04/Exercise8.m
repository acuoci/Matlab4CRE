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
% Adiabatic CSTR                                                          %
% Two parallel first order reactions                                      %
% Direct solution of a Non Linear System of algebraic equations           %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all;
clear all;

global A1;           % [1/s]
global Ea1;          % [cal/mol]
global A2;           % [1/s]
global Ea2;          % [cal/mol]
global CAin;         % [kmol/m3]
global CBin;         % [kmol/m3]
global CCin;         % [kmol/m3]
global tau;          % [s]
global Tin;          % [K]
global T0;           % [K]
global HA0;          % [kcal/kmol]
global HB0;          % [kcal/kmol]
global HC0;          % [kcal/kmol]
global CpA;          % [kcal/kmol/K]
global CpB;          % [kcal/kmol/K]
global CpC;          % [kcal/kmol/K]

% ------------------------------------------------------------------------%
% User data                                                               %
% ------------------------------------------------------------------------%
A1 = 3e14;              % frequency factor [1/s]
Ea1 = 20000;            % activation energy [cal/mol]
A2 = 2e13;              % frequency factor [1/s]
Ea2 = 18000;            % activation temperature [cal/mol]
Tin = 20+273.15;        % reactor temperature [K]
tau = 1*60;             % residence time [s]
CAin = 55;              % inlet concentration of A [kmol/m3]
CBin = 0;               % inlet concentration of B [kmol/m3]
CCin = 0;               % inlet concentration of C [kmol/m3]
T0 = 20+273.15;         % reference temperature [K]
HA0 =  6000;            % enthalpy of species A at T0 [kcal/kmol]
HB0 =  5802;            % enthalpy of species B at T0 [kcal/kmol]
HC0 =  5620;            % enthalpy of species C at T0 [kcal/kmol]
CpA = 5;                % specific heat of A [kcal/kmol/K]
CpB = 6;                % specific heat of B [kcal/kmol/K]
CpC = 4;                % specific heat of C [kcal/kmol/K]

% ------------------------------------------------------------------------%
% NLS solution                                                            %
% ------------------------------------------------------------------------%
YFirstGuess = [CAin CBin CCin Tin+10]';
options = optimset('Display', 'off');
Y = fsolve(@CSTR, YFirstGuess, options);
CA = Y(1);
CB = Y(2);
CC = Y(3);
T  = Y(4)-273.15;

% Plot results
fprintf('A) CAout: %f [kmol/m3]\n', CA);
fprintf('B) CBout: %f [kmol/m3]\n', CB);
fprintf('C) CCout: %f [kmol/m3]\n', CC);
fprintf('T)  Tout: %f [C]\n', T);

% ------------------------------------------------------------------------%
% NLS definition                                                           %
% ------------------------------------------------------------------------%
function F = CSTR(Y)

    global A1;           % [1/s]
    global Ea1;          % [cal/mol]
    global A2;           % [1/s]
    global Ea2;          % [cal/mol]
    global CAin;         % [kmol/m3]
    global CBin;         % [kmol/m3]
    global CCin;         % [kmol/m3]
    global tau;          % [s]
    global Tin;          % [K]
    global T0;           % [K]
    global HA0;          % [kcal/kmol]
    global HB0;          % [kcal/kmol]
    global HC0;          % [kcal/kmol]
    global CpA;          % [kcal/kmol/K]
    global CpB;          % [kcal/kmol/K]
    global CpC;          % [kcal/kmol/K]
    
    % Recover unknowns 
    CA = Y(1);                  % [kmol/m3]
    CB = Y(2);                  % [kmol/m3]
    CC = Y(3);                  % [kmol/m3]
    T  = Y(4);                  % [kmol/m3]
    
    % Reaction rates
    k1 = A1*exp(-Ea1/1.987/T);  % [1/s]
    k2 = A2*exp(-Ea2/1.987/T);  % [1/s]
    r1 = k1*CA;                 % [kmol/m3/s]
    r2 = k2*CA;                 % [kmol/m3/s]
    
    % Formation rates
    RA = -r1-r2;                % [kmol/m3/s]
    RB = r1;                    % [kmol/m3/s]
    RC = r2;                    % [kmol/m3/s]
    
    % Enthalpies
    HAin = HA0+CpA*(Tin-T0);
    HBin = HB0+CpB*(Tin-T0);
    HCin = HC0+CpC*(Tin-T0);
    HA   = HA0+CpA*(T-T0);
    HB   = HB0+CpB*(T-T0);
    HC   = HC0+CpC*(T-T0);
    
    % Non linear equations
    F(1) = (CAin-CA)/tau + RA;
    F(2) = (CBin-CB)/tau + RB;
    F(3) = (CCin-CC)/tau + RC;
    F(4) = HAin*CAin+HBin*CBin+HCin*CCin - (HA*CA+HB*CB+HC*CC);
    
end