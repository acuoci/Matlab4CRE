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
%   Copyright(C) 2018 Alberto Cuoci                                       |
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
close all;
clear all;

global pTot;        % total pressure [Pa]
global Ctot;        % total concentration liquid phase [kmol/m3]
global FtotG;       % total molar flow rate (gas) [kmol/s]
global FtotL;       % total molar flow rate (liquid) [kmol/s]
global b;           % stoichiometric coefficient species B 
global c;           % stoichiometric coefficient species C

b = 2;          % stoichiometric coefficient species B
c = 1;          % stoichiometric coefficient species C
T = 303;        % temperature [K]
pTot = 101325;  % total pressure (gas) [Pa]
pAin = 5000;    % inlet partial pressure of species A [Pa]
CBin = 0.001;   % inlet concentration of species B [kmol/m3]
CCin = 0.0;     % inlet concentration of species C [kmol/m3]
QL = 20/1000;   % liquid volumetric flow rate [m3/s]
Ctot = 0.301;   % total concentration liquid phase [kmol/m3]
QG = 0.8/1000;  % gas volumetric flow rate [m3/s]
L = 5;          % total length [m]
D = 0.30;       % diameter [m]

% ------------------------------------------------------------------------%
% Co-current Tower Reactor                                                %
% ------------------------------------------------------------------------%
FtotG = (pTot/8314/T)*QG;   % total molar flow rate (gas) [kmol/s]
FtotL = Ctot*QL;            % total molar flow rate (liquid) [kmol/s]
A = pi/4*D^2;               % cross section area [m2]
Vtot = A*L;                 % total volume [m3]

% ODE solution
Yin = [pAin CBin CCin]';
[V, Y] = ode45(@TowerCoCurrent, [0 Vtot], Yin);
z = V/A;        % axial coordinate [m]
pA = Y(:,1);    % partial pressure of A [Pa]
CB = Y(:,2);    % concentration of B [kmol/m3]
CC = Y(:,3);    % concentration of C [kmol/m3]

% ------------------------------------------------------------------------%
% Post-processing
% ------------------------------------------------------------------------%

fprintf('A conversion: %f \n', 1-pA(end)/pAin);

figure;
plot(z,pA);
xlabel('axial coordinate [m]'); ylabel('partial pressure of A [Pa]'); 
title('gaseous phase');

figure;
plot(z,CB*1000);
xlabel('axial coordinate [m]'); ylabel('concentration of B [mol/m3]'); 
title('liquid phase');

figure;
plot(z,CC*1000);
xlabel('axial coordinate [m]'); ylabel('concentration of C [mol/m3]'); 
title('liquid phase');

% Reconstructing the overall rate of change and the enhancing factor
for i=1:size(pA,1)
    [r(i),E(i)] = OverallRateOfChange(pA(i), CB(i));
end

figure;
plot(z,r*3600*1000);
xlabel('axial coordinate [m]'); ylabel('overall rate of change [mol/m3/hr]'); 
title('overall rate of change ');

figure;
plot(z,E);
xlabel('axial coordinate [m]'); ylabel('enhancing factor'); 
title('enhancing factor ');


% ------------------------------------------------------------------------%
% ODE system (defined locally, requires at least MATLAB-2016b)            %
% ------------------------------------------------------------------------%
function dY = TowerCoCurrent(V,Y)

    global pTot;    % total pressure [Pa]
    global Ctot;    % total concentration [kmol/m3]
    global FtotG;   % total molar flow rate (gas) [kmol/s]
    global FtotL;   % total molar flow rate (liquid) [kmol/s]
    global b;       % stoichiometric coefficient B
    global c;       % stoichiometric coefficient C

    pA = Y(1);      % partial pressure of A in gaseous phase [Pa]
    CB = Y(2);      % concentration of B in liquid phase [kmol/m3]
    
    [r] = OverallRateOfChange(pA, CB);      % [kmol/m3/s]
        
    dpA = -pTot/FtotG*r;                    % [Pa/m3]
    dCB = -Ctot/FtotL*b*r;                  % [kmol/m3/m3]
    dCC =  Ctot/FtotL*c*r;                  % [kmol/m3/m3]

    dY = [dpA dCB dCC]';

end
