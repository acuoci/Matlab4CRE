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
% Cocurrent bubble tower                                                  %
% Third order reaction: A+2B->C                                           %
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
data.kIII=100e6/1e-6/3600.;      % kinetic constant (liquid volume basis) [m6/kmol^2/s]


%% Input data
T = 303;        % temperature [K]
pTot = 101325;  % total pressure (gas) [Pa]
pAin = 5000;    % inlet partial pressure of species A [Pa]
CBin = 0.001;   % inlet concentration of species B [kmol/m3]
CCin = 0.0;     % inlet concentration of species C [kmol/m3]
QL = 10/1000;   % liquid volumetric flow rate [m3/s]
Ctot = 0.301;   % total concentration liquid phase [kmol/m3]
QG = 0.4/1000;  % gas volumetric flow rate [m3/s]
L = 5;          % total length [m]
D = 0.40;       % diameter [m]


%% Co-current Tower Reactor
FtotG = (pTot/8314/T)*QG;   % total molar flow rate (gas) [kmol/s]
FtotL = Ctot*QL;            % total molar flow rate (liquid) [kmol/s]
A = pi/4*D^2;               % cross section area [m2]
Vtot = A*L;                 % total volume [m3]


%% ODE solution
Yin = [pAin CBin CCin]';
[V, Y] = ode45(@TowerCoCurrent, [0 Vtot], Yin, [], pTot, Ctot, FtotG, FtotL, data);
z = V/A;        % axial coordinate [m]
pA = Y(:,1);    % partial pressure of A [Pa]
CB = Y(:,2);    % concentration of B [kmol/m3]
CC = Y(:,3);    % concentration of C [kmol/m3]


%% Post-processing

fprintf('A conversion: %f \n', 1-pA(end)/pAin);

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('partial pressure of A [Pa]'); 
plot(z,pA);
yyaxis right; ylabel('concentration of B [mol/m3]'); 
plot(z,CB*1000);

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('concentration of B [mol/m3]'); 
plot(z,CB*1000);
yyaxis right; ylabel('concentration of C [mol/m3]'); 
plot(z,CC*1000);


%% Reconstructing the overall rate of change and the enhancing factor
for i=1:size(pA,1)
    [r(i),E(i)] = OverallRateOfChange(pA(i), CB(i), data);
end

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('overall rate of change [mol/m3/hr]'); 
plot(z,r*3600*1000);
yyaxis right; ylabel('enhancing factor'); 
plot(z,E);


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = TowerCoCurrent(V,Y, pTot, Ctot, FtotG, FtotL, data)

    pA = Y(1);      % partial pressure of A in gaseous phase [Pa]
    CB = Y(2);      % concentration of B in liquid phase [kmol/m3]
    
    [r] = OverallRateOfChange(pA, CB, data);    % [kmol/m3/s]
        
    dpA = -pTot/FtotG*r;                        % [Pa/m3]
    dCB = -Ctot/FtotL*data.b*r;                 % [kmol/m3/m3]
    dCC =  Ctot/FtotL*data.c*r;                 % [kmol/m3/m3]

    dY = [dpA dCB dCC]';

end
