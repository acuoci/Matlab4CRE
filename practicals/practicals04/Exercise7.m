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
% Isothermal CSTR                                                         %
% Two parallel first order reactions                                      %
% Analytical solution available                                           %
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
global T;            % [K]

% ------------------------------------------------------------------------%
% User data                                                               %
% ------------------------------------------------------------------------%
A1 = 3e14;              % frequency factor [1/s]
Ea1 = 20000;            % activation energy [cal/mol]
A2 = 2e13;              % frequency factor [1/s]
Ea2 = 18000;            % activation temperature [cal/mol]
T = 80+273.15;          % reactor temperature [K]
tau = 1*60;             % residence time [s]
CAin = 55;              % inlet concentration of A [kmol/m3]
CBin = 0;               % inlet concentration of A [kmol/m3]
CCin = 0;               % inlet concentration of A [kmol/m3]

% ------------------------------------------------------------------------%
% NLS solution                                                            %
% ------------------------------------------------------------------------%
YFirstGuess = [CAin CBin CCin]';
options = optimset('Display', 'off');
Y = fsolve(@CSTR, YFirstGuess, options);
CAnum = Y(1);
CBnum = Y(2);
CCnum = Y(3);

% Analytical Solution
k1 = A1*exp(-Ea1/1.987/T);  % [1/s]
k2 = A2*exp(-Ea2/1.987/T);  % [1/s]
ktot = k1+k2;
CAexact = CAin/(1+ktot*tau);
CBexact = CBin +k1*tau*CAexact;
CCexact = CCin +k2*tau*CAexact;

% Comparison
fprintf('A) Numerical: %f - Exact: %f [kmol/m3]\n', CAnum, CAexact);
fprintf('B) Numerical: %f - Exact: %f [kmol/m3]\n', CBnum, CBexact);
fprintf('C) Numerical: %f - Exact: %f [kmol/m3]\n', CCnum, CCexact);

% ------------------------------------------------------------------------%
% NLS definition                                                          %
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
    global T;            % [K]
    
    % Recover concentrations 
    CA = Y(1);                  % [kmol/m3]
    CB = Y(2);                  % [kmol/m3]
    CC = Y(3);                  % [kmol/m3]
    
    % Reaction rates
    k1 = A1*exp(-Ea1/1.987/T);  % [1/s]
    k2 = A2*exp(-Ea2/1.987/T);  % [1/s]
    r1 = k1*CA;                 % [kmol/m3/s]
    r2 = k2*CA;                 % [kmol/m3/s]
    
    % Formation rates
    RA = -r1-r2;                % [kmol/m3/s]
    RB = r1;                    % [kmol/m3/s]
    RC = r2;                    % [kmol/m3/s]
    
    % Non linear equations
    F(1) = (CAin-CA)/tau + RA;
    F(2) = (CBin-CB)/tau + RB;
    F(3) = (CCin-CC)/tau + RC;
    
end