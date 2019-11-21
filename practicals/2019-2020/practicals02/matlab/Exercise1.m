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
% Non isothermal plug flow reactor with heat exchange                     %
% Constant density, first order reaction                                  %
% Exothermic reaction, constant temperature external fluid                %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all; clear variables;

global A;       % [1/s]
global Ea;      % [cal/mol]
global P;       % [Pa]
global S;       % [m2]
global Te;      % [K]
global U;       % [cal/m2/h/K]
global D;       % [m]
global deltaHR; % [cal/mol]
global Cp;      % [cal/mol/K]


%% User data
A = 2e8;            % frequency factor [1/s]
Ea = 24000;         % activation energy [cal/mol]
deltaHR = -5000;    % reaction heat [cal/mol]
Cp = 30;            % cp [cal/mol/K]

Te = 300 + 273.15;  % external temperature [K]
P = 3*101325;       % pressure [Pa]

FAin = 50;          % inlet flow rate [kmol/h]
FBin = 0;           % inlet flow rate [kmol/h]
Tin = 300 + 273.15; % inlet temperature [K]

U = 50;             % global heat exchange coefficient [kcal/m2/h/K]
D = 8e-2;           % internal diameter [m]
L = 150;            % reactor length [m]


%% Calculations

% Cross section area
S = pi/4*D^2;   % cross section area [m2]

%ODE solution
Yin = [FAin FBin Tin]';
[z, Y] = ode45(@PFR, [0 L], Yin);


%% Post processing
FA = Y(:,1);                % molar flow rate [kmol/h]
FB = Y(:,2);                % molar flow rate [kmol/h]
T  = Y(:,3)-273.15;         % temperature [C]
X = (FAin - FA) / FAin;     % conversion

%Plotting results
figure;
plot(z,X);
xlabel('axial coordinate [m]'); ylabel('conversion'); title('conversion');

figure;
plot(z,T);
xlabel('axial coordinate [m]'); ylabel('T [C]'); title('temperature');

figure;
plot(z,FA, z,FB);
xlabel('axial coordinate [m]'); ylabel('molar flow rates [kmol/h]'); 
title('molar flow rates'); legend('A', 'B');

[maxT, iMaxT] = max(T);
fprintf('Exit: Conversion: %f Temperature: %f C\n', X(end), T(end));
fprintf('Max temperature: %f C @ z=%f\n', maxT, z(iMaxT));


%% Parametric analysis
figure;
hold on;
xlabel('axial coordinate [m]'); ylabel('T [C]'); title('temperature');

for i=1:6
    
    U = (i-1)*20;
    [z, Y] = ode45(@PFR, [0 L], Yin);
    T = Y(:,3)-273.15;
    plot(z,T);
    legendInfo{i} = ['U=' num2str(U)];
    
end
refline(0,380);
legendInfo{7}='Reference'; legend(legendInfo, 'Location','northwest');


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = PFR(z,Y)

    global A;       % [1/s]
    global Ea;      % [cal/mol]
    global P;       % [Pa]
    global S;       % [m2]
    global Te;      % [K]
    global U;       % [cal/m2/h/K]
    global D;       % [m]
    global deltaHR; % [cal/mol]
    global Cp;      % [cal/mol/K]

    FA = Y(1);      % molar flow rate of A [kmol/h]
    FB = Y(2);      % molar flow rate of B [kmol/h]
    T  = Y(3);      % molar flow rate of C [kmol/h]
    
    Ctot = P/8314/T;                % [kmol/m3]
    k = A*exp(-Ea/1.987/T)*3600;    % [1/h]
    
    Ftot = FA+FB;                   % [kmol/h]
    CA = Ctot*FA/Ftot;              % [kmol/m3]
    r = k*CA;                       % [kmol/m3/h]
    
    dFA = -S*r;                                     % [kmol/m/h]
    dFB =  S*r;                                     % [kmol/m/h]
    dT  = -S*( U*(T-Te)*4/D + r*deltaHR)/(Ftot*Cp); % [K/h]

    dY = [dFA dFB dT]';

end
