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
% Adiabatic PFR reactor (constant density)                                %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all;  clear variables;


%% Input data
kfRef=31.1;  % [1/h]
Tfref=360;   % [K]
E=65700;     % [J/mol]		
KeqRef=3.03; % [-]	
TeqRef=333;  % [K]		
DH0=-6900;   % [J/mol]
T0=300;      % [K]
CpA=131;     % [J/mol/K]
CpB=171;     % [J/mol/K]
CpI=161;     % [J/mol/K]	
CAin=9.3;    % [kmol/m3]
FtotIn=163;  % [kmol/h]
xAin=0.9;    % [-]	
xBin=0;      % [-]
xIin=0.1;    % [-]	
Tin=330;     % [K]		
XTarget=0.6; %[-]


%% Calculations

% Thermodynamic properties
tetaA=xAin/xAin;
tetaB=xBin/xAin;
tetaI=xIin/xAin;
Cpin=tetaA*CpA+tetaB*CpB+tetaI*CpI;
DCp=CpB-CpA;
DHStar=DH0+DCp*(TeqRef-T0);
DHin=DH0+DCp*(Tin-T0);
FAin = FtotIn*xAin;    % [kmol/h]

%% ODE solution
Y0 = [0]';
opt_ode = odeset('RelTol',1e-7,'AbsTol',1e-12);
[V, Y] = ode15s(@adiabaticPFR, [0:0.01:2.5], Y0, opt_ode, ...
                kfRef,Tfref,E,KeqRef,TeqRef,DCp,DHStar,DH0,T0,CAin,FAin,Cpin,Tin);

% Post processing
X = Y(:,1);                % [kmol/m3]

% Conversion
figure;
plot(V,X,'-o');
xlabel('volume [m3]');
ylabel('conversion');
title('Conversion');

% Find volume
iTarget = find(X>=XTarget,1);
VTarget = V(iTarget);               % [m3]
fprintf('Target: V(X=%f): %f \n', X(iTarget), VTarget);


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = adiabaticPFR(V,Y, kfRef,Tfref,E,KeqRef,TeqRef,...
                           DCp,DHStar,DH0,T0,CAin,FAin,Cpin,Tin)
    
    % Recover variables
    X = Y(1);      % [-]

    % Reconstruct temperature
    DHin=DH0+DCp*(Tin-T0);
    T = Tin + (-DHin*X)/(Cpin+DCp*X);
    
    % Kinetic data
    kf  = kfRef*exp(-E/8.314*(1/T-1/Tfref));    % [1/h]
    kEq = KeqRef*exp(-(DHStar/8.314-DCp/8.314*TeqRef)*(1/T-1/TeqRef) + ...
          DCp/8.314*log(T/TeqRef));
    r  = kf*CAin*(1-(1+1/kEq)*X);               % [kmol/m3/h]
    
    % Equations
    dX = r/FAin;                 % [1/m3]
    dY = [dX]';

end
