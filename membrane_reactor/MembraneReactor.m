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
% Comparison PFR vs Membrane Reactor: 
%
%              A <-> B + C
%
% For details see: Fogler, Chapter 4 (Web modules: membrane reactors)
%
%-------------------------------------------------------------------------%
close all;
clear all;

% ------------------------------------------------------------------------%
% Global data
% ------------------------------------------------------------------------%
global k;
global KC;
global T0;
global P0;
global km;

% ------------------------------------------------------------------------%
% Definition of functions corresponding to ODE systems to be solved
% ------------------------------------------------------------------------%
myfuns = MembraneReactorFunctions;

% ------------------------------------------------------------------------%
% Data
% ------------------------------------------------------------------------%
P0   = 6*101325; % pressure [Pa]
T0   = 373;      % temperature [K]
km   = 0.2;      % mass transfer coefficient [1/min]
FA0  = 15;       % inlet flow rate of species A [mol/min]
FB0  =  0;       % inlet flow rate of species B [mol/min]
FC0  =  0;       % inlet flow rate of species C [mol/min]
dV = 1000;       % final volume [dm3]

% ------------------------------------------------------------------------%
% Kinetic and thermodynamic data
% ------------------------------------------------------------------------%
k = 0.7;   % forward kinetic constant [1/min]                  
KC = 0.05; % equilibrium constant[mol/dm3]  

% ------------------------------------------------------------------------%
% Initial conditions
% ------------------------------------------------------------------------%
y0 = [FA0, FB0, FC0];

%*************************************************************************%
% Part 1: pure plug flow reactor
%*************************************************************************%

% ------------------------------------------------------------------------%
% ODE solution
% ------------------------------------------------------------------------%
[V,y] = ode15s(@myfuns.plugflow,[0 dV], y0);

% ------------------------------------------------------------------------%
% Plotting solution
% ------------------------------------------------------------------------%
figure;
subplot(2,3,1);
hold all;
plot(V,y(:,1), 'LineWidth', 2);
plot(V,y(:,2), 'LineWidth', 2);
plot(V,y(:,3), 'LineWidth', 2, 'LineStyle', '--');
title ('Pure plug flow reactor');
legend('A','B','C');
xlabel('volume [dm3]');
ylabel('flow rates [mol/min]');

subplot(2,3,2);
hold all;
Ctot = P0/8314/T0;
Ftot = y(:,1)+y(:,2)+y(:,3);
plot(V,Ctot*y(:,1)./Ftot, 'LineWidth', 2);
plot(V,Ctot*y(:,2)./Ftot, 'LineWidth', 2);
plot(V,Ctot*y(:,3)./Ftot, 'LineWidth', 2, 'LineStyle', '--');
title ('Pure plug flow reactor');
legend('A','B','C');
xlabel('volume [dm3]');
ylabel('concentrations [mol/dm3]');

subplot(2,3,3);
hold all;
[hAx,hLine1,hLine2] = plotyy(V,Ftot, V, (FA0-y(:,1))/FA0);
title ('Pure plug flow reactor');
legend('Total molar flow rate', 'Conversion');
xlabel('volume [dm3]');
set(hLine1, 'LineWidth', 2);
set(hLine2, 'LineWidth', 2);
ylabel(hAx(1),'total molar flow rate [mol/min]') % left y-axis
ylabel(hAx(2),'conversion [-]')                  % right y-axis

%*************************************************************************%
% Part 2: membrane reactor (only B component)
%*************************************************************************%

% ------------------------------------------------------------------------%
% ODE solution
% ------------------------------------------------------------------------%
[V,y] = ode15s(@myfuns.membranereactor_only_B,[0 dV], y0);

% ------------------------------------------------------------------------%
% Plotting solution
% ------------------------------------------------------------------------%
subplot(2,3,4);
hold all;
plot(V,y(:,1), 'LineWidth', 2);
plot(V,y(:,2), 'LineWidth', 2);
plot(V,y(:,3), 'LineWidth', 2, 'LineStyle', '--');
title ('Membrane reactor (only B)');
legend('A','B','C');
xlabel('volume [dm3]');
ylabel('flow rates [mol/min]');

subplot(2,3,5);
hold all;
Ctot = P0/8314/T0;
Ftot = y(:,1)+y(:,2)+y(:,3);
plot(V,Ctot*y(:,1)./Ftot, 'LineWidth', 2);
plot(V,Ctot*y(:,2)./Ftot, 'LineWidth', 2);
plot(V,Ctot*y(:,3)./Ftot, 'LineWidth', 2, 'LineStyle', '--');
title ('Membrane reactor (only B)');
legend('A','B','C');
xlabel('volume [dm3]');
ylabel('concentrations [mol/dm3]');

subplot(2,3,6);
hold all;
[hAx,hLine1,hLine2] = plotyy(V,Ftot, V, (FA0-y(:,1))/FA0);
title ('Pure plug flow reactor');
legend('Total molar flow rate', 'Conversion');
title ('Membrane reactor (only B)');
xlabel('volume [dm3]');
set(hLine1, 'LineWidth', 2);
set(hLine2, 'LineWidth', 2);
ylabel(hAx(1),'total molar flow rate [mol/min]') % left y-axis
ylabel(hAx(2),'conversion [-]')                  % right y-axis


