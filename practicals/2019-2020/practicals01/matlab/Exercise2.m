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
% Batch reactor (constant density)                                        %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all;  clear variables;


%% Input data
T = 100+273.15; % [K]
V = 1;          % [m3]
TauD = 1;       % [h]
rho = 800;      % [kg/m3]
MW = 35;        % [kg/kmol]
XTarget = 0.98; % [-]
PriceB = 15;    % [$/kg]
C1 = 100;       % [$/h]
C2 = 25;        % [$/h]
C3 = 8000;      % [$/h]

% Kinetic constants
k1 = 8e4*exp(-8000/1.987/T);   % [1/h]
k2 = 1e5*exp(-10000/1.987/T);  % [1/h]

% Initial conditions
Ctot = rho/MW;          % [kmol/m3]
CA0  = Ctot;            % [kmol/m3]
CB0  = 0;               % [kmol/m3]
CC0  = 0;               % [kmol/m3]
tau  = 3;               % [h] (large value)


%% ODE solution
C0 = [CA0 CB0 CC0]';
[t, C] = ode45(@isothermalBatch, [0 tau], C0, [], k1,k2);


%% Post processing
CA = C(:,1);        % [kmol/m3]
CB = C(:,2);        % [kmol/m3]
CC = C(:,3);        % [kmol/m3]
X = (CA0-CA)/CA0;   % [-]

% Concentrations
figure;
plot(t,CA,'-o', t,CB,'-o', t,CC,'-o');
legend('CA', 'CB', 'CC');
xlabel('time [h]');
ylabel('C [kmol/m3]');
title('Concentrations');

% Conversion
figure;
plot(t,X,'-o');
xlabel('time [h]');
ylabel('X');
title('Conversion');

% Find residence time (target conditions)
iTarget = find(X>=XTarget,1);
tauTarget = t(iTarget);             % [h]
NTarget = 24/(tauTarget+TauD);      % [#/day]
PBTarget = CB(iTarget)*V*NTarget;   % [kmol/day]

fprintf('Target: Tau(X=%f): %f \n', XTarget, tauTarget);
fprintf('Target: N (per day): %f \n', NTarget);
fprintf('Target: PB: %f [kmol/day]\n', PBTarget);

% Optimal conditions
Incomes = PriceB*(CB*MW)*V*24./(t+TauD);    % [$]
Costs = (C1*TauD+C2*t+C3)*24./(t+TauD);     % [$]
Margin = Incomes-Costs;                     % [$]

% Economic analysis
figure;
plot(t,Incomes,'-o', t,Costs,'-o', t,Margin,'-o');
legend('I', 'C', 'M');
xlabel('time [h]');
ylabel('$');
title('Economic Analysis');

[maxMargin, iMaxMargin] = max(Margin);
tauMaxMargin = t(iMaxMargin);
fprintf('Optimization: Max Margin($): %f @ Tau(s)=%f \n', ...
         maxMargin, tauMaxMargin);


%% ODE system (defined locally, requires at least MATLAB-2016b
function dC = isothermalBatch(t,C, k1,k2)

    CA = C(1);      % [kmol/m3]
    CB = C(2);      % [kmol/m3]
    CC = C(3);      % [kmol/m3]

    r1 = k1*CA;     % [kmol/m3/s]
    r2 = k2*CA;     % [kmol/m3/s]

    dCA = -r1-r2;   % [kmol/m3/h]
    dCB = r1;       % [kmol/m3/s]
    dCC = r2;       % [kmol/m3/s]

    dC = [dCA dCB dCC]';

end
