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

%-------------------------------------------------------------------------%
% Fogler, Exercise 3-6:  Calculating the equilibrium conversion
%
% Reaction: A <-> 2B 
% r = k*(CA - CB^2/KC)
%-------------------------------------------------------------------------%
clear all;
	
%-------------------------------------------------------------------------%
% Data
%-------------------------------------------------------------------------%
eps = 1;		% changes of moles
P = 202.6e3;	% pressure [Pa]
T = 340;		% temperature [K]
KC = 0.1;		% equilibrium constant [mol/dm3]
k1 = 0.5;		% forward kinetic constant [1/s]
CA0 = P/8314/T; % initial concentration [kmol/m3]=[mol/dm3]


%-------------------------------------------------------------------------%
% 1. Batch Reactor (BR)
%-------------------------------------------------------------------------%
equilibriumXBatch = @(x) 4*CA0*x^2+x*KC-KC;  % parametrized function
fun = @(x) equilibriumXBatch(x);    		 % function of x alone
xBatch = fzero(fun,0.2);

XBatch = linspace(0,xBatch);
batch_minus_RA = k1*CA0*((1-XBatch)-4*CA0*XBatch.^2/KC);
batch_CA = CA0*(1-XBatch);
batch_CB = 2*CA0*XBatch;

%-------------------------------------------------------------------------%
% 2. Flow Reactors
%-------------------------------------------------------------------------%
equilibriumXFlow = @(x) 4*CA0*x^2-KC*(1-x)*(1+eps*x);  	% parametrized function
fun = @(x) equilibriumXFlow(x);    						% function of x alone
xFlow = fzero(fun,0.2);

XFlow = linspace(0,xFlow);
flow_minus_RA = k1*CA0./(1+eps*XFlow).*((1-XFlow)-4*CA0*XFlow.^2/KC./(1+eps*XFlow));
flow_CA = CA0*(1-XFlow)./(1+eps*XFlow);
flow_CB = 2*CA0*XFlow./(1+eps*XFlow);

%-------------------------------------------------------------------------%
% Plots
%-------------------------------------------------------------------------%
figure; % create new figure

subplot(2,2,1);
hold all;
plot(XBatch,1./batch_minus_RA, 'LineWidth',2);
plot(XFlow,1./flow_minus_RA,'LineWidth',2);
title('Levenspiel Plot');
xlim([0 0.5]);
ylim([0 1000]);
xlabel('conversion');
ylabel('-1/R_A [dm^3.min/mol]');
legend('Batch','Flow');

subplot(2,2,2);
hold all;
plot(XBatch,batch_minus_RA,'LineWidth',2);
plot(XFlow,flow_minus_RA,'LineWidth',2);
title('Reaction rate plot');
xlim([0 0.5]);
%ylim([0 1000]);
xlabel('conversion');
ylabel('r [mol/dm^3/min]');
legend('Batch','Flow');

subplot(2,2,3);
hold all;
plot(XBatch,batch_CA,'LineWidth',2);
plot(XFlow,flow_CA,'LineWidth',2);
title('Concentration species A');
xlim([0 0.5]);
xlabel('conversion');
ylabel('CA [mol/dm^3]');
legend('Batch','Flow');

subplot(2,2,4);
hold all;
plot(XBatch,batch_CB,'LineWidth',2);
plot(XFlow,flow_CB,'LineWidth',2);
title('Concentration species B');
xlim([0 0.5]);
xlabel('conversion');
ylabel('CB [mol/dm^3]');
legend('Batch','Flow');
