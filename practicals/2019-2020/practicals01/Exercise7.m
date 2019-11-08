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
%   Copyright(C) 2019 Alberto Cuoci                                       |
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
% Isothermal, batch reactor (constant density)                            %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all;	clear variables;

%% Global variables
global k1;
global k2;


%% Data
A1=1.50E+04;	% [1/h]
E1=9000;        % [cal/mol]
A2=6.00E+06;	% [1/h]
E2=19000;       % [cal/mol]
		
T=500;	% [K]
V=0.5;	% [m3]
NA0=20;	% [kmol]
TauD=1;	% [h]


%% Preliminary calculations

% Kinetic constants
k1=A1*exp(-E1/1.987/T);  %[1/h]
k2=A2*exp(-E2/1.987/T);  %[1/h]

% Initial conditions
CA0 = NA0/V;    % [kmol]
CB0 = 0.;       % [kmol]
CC0 = 0.;       % [kmol]
tau = 3;        % [h]


%% ODE solution
C0 = [CA0 CB0 CC0]';
[t, C] = ode45(@isothermalBatch, [0 tau], C0);


%% Post processing
CA = C(:,1);                % [kmol/m3]
CB = C(:,2);                % [kmol/m3]
CC = C(:,3);                % [kmol/m3]
PB = CB*V*24./(TauD+t);     % [kmol/day]

% Concentrations
figure;
plot(t,CA,'-o', t,CB,'-o', t,CC,'-o');
legend('CA', 'CB', 'CC');
xlabel('time [h]');
ylabel('C [kmol/m3]');
title('Concentrations');

[maxCB, iMaxCB] = max(CB);
tauMaxCB = t(iMaxCB);
fprintf('Optimization: Max CB($): %f @ Tau(s)=%f \n', ...
         maxCB, tauMaxCB);

%Production
figure;
plot(t,PB,'-o');
xlabel('time [h]');
ylabel('PB [kmol/day]');
title('Production');

[maxPB, iMaxPB] = max(PB);
tauMaxPB = t(iMaxPB);
fprintf('Optimization: Max PB($): %f @ Tau(s)=%f \n', ...
         maxPB, tauMaxPB);



%% ODE system (defined locally, requires at least MATLAB-2016b)
function dC = isothermalBatch(t,C)

    global k1;      % [1/h]
    global k2;      % [1/h]

    CA = C(1);      % [kmol/m3]
    CB = C(2);      % [kmol/m3]

    r1 = k1*CA;     % [kmol/m3/h]
    r2 = k2*CB;     % [kmol/m3/h]

    dCA = -r1;      % [kmol/m3/h]
    dCB = r1-r2;    % [kmol/m3/h]
    dCC = r2;       % [kmol/m3/h]

    dC = [dCA dCB dCC]';

end
