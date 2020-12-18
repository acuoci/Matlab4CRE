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
%   Copyright(C) 2020 Alberto Cuoci                                       |
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
% Calculation of overall rate of change (gas/liquid)                      %
%                                                                         %
%-------------------------------------------------------------------------%
close all; clear variables;

%% Data struct
data.b=2;                        % stoichiometric coefficient species B
data.c=1;                        % stoichiometric coefficient species C
data.GammaA=1e-6/3600;           % diffusion coefficient [m2/s]
data.GammaB=1e-6/3600;           % diffusion coefficient [m2/s]
data.fL=0.98;                    % liquid fraction [-]
data.HA=1e5*1e3;                 % Henry's constant [Pa.m3/kmol]
data.a=20;                       % interfacial area per unti ov volume [m2/m3]
data.KLa=20/3600.;               % mass transfer coefficient (liquid side) [1/s]
data.KGa=0.01*(1e-3/3600);       % mass transfer coefficient (gas side) [kmol/m3/s/Pa]
data.kl=100e6/1e-6/3600.;        % kinetic constant (liquid volume basis) [m6/kmol^2/s]

%% Input data
pTot = 101325;  % total pressure (gas) [Pa]
Ctot = 0.301;   % total concentration liquid phase [kmol/m3]
pAin = 5000;    % inlet partial pressure of species A [Pa]
CBin = 0.001;   % inlet concentration of species B [kmol/m3]
QL = 20/1000;   % liquid volumetric flow rate [m3/s]
QG = 0.8/1000;  % gas volumetric flow rate [m3/s]
T = 303;        % temperature [K]

%% Pre-processing
FtotG = (pTot/8314/T)*QG;   % total molar flow rate (gas) [kmol/s]
FtotL = Ctot*QL;            % total molar flow rate (liquid) [kmol/s]

%% Partial pressure of species A [Pa]
pA = 100:100:pAin;  

%% Operating line: Concentration of species B [kmol/m3]
CB = CBin + FtotG/FtotL*data.b*Ctot/pTot*(pA-pAin);

%% Rate of change [kmol/m3/s]
for i=1:length(pA)
    [rIIII(i),E(i),MH(i)] = OverallRateOfChange(pA(i), CB(i), data);
end

%% Levenspiel's plots
figure;
plot(pA,1./(rIIII*3600*1000));
xlabel('partial pressure of A [Pa]'); 
ylabel('1/rIIII [m3.hr/mol]'); 
title('Levenspiel Plot');

figure;
plot(CB*1000,1./(rIIII*3600*1000));
xlabel('concentration of B [mol/m3]'); 
ylabel('1/rIIII [m3.hr/mol]'); 
title('Levenspiel Plot');

%% Operating line
figure;
plot(CB*1000,pA);
xlabel('concentration of B [mol/m3]'); 
ylabel('partial pressure of A [Pa]'); 
title('Operating line');

%% Enhancement factor
figure;
plot(pA,E);
xlabel('partial pressure of A [Pa]'); 
ylabel('Enhancement factor'); 
