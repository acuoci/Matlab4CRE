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

clear all;

% ------------------------------------------------------------------------%
% Global variables 
% ------------------------------------------------------------------------%
global kinetics;            % kinetic mechanism
global T0;                  % inlet temperature [K]
global P0;                  % inlet pressure [Pa]
global omega0;              % inlet mass fractions [-]
global v0;                  % inlet velocity [m/s]
global D;                   % internal diameter
global gx;                  % gravity [m/s2]
global U;
global Te;

% ------------------------------------------------------------------------%
% 1. Load the kinetic mechanism 
% ------------------------------------------------------------------------%
kinetics = GasMechanism_AceticAnhydride;

% ------------------------------------------------------------------------%
% 2. Define input data
% ------------------------------------------------------------------------%
T0 = 1000.;                     % inlet temperature [K]
P0 = 101325;                    % inlet pressure [Pa]

omega0 = zeros(1,kinetics.ns);  % inlet mass fractions [-]
omega0(1) = 1.;

D = 2.66e-2;                    % internal diameter [m]
v0 = 10;                        % inlet velocity [m/s]
x0 = 0;                         % initial axial coordinate [m]
xF = 10;                        % final axial coordinate [m]

U  = 30;                        % global heat exchange coefficient [W/m2/K]
Te = 1100;                      % external temperature [K]

gx = 9.81;                      % gravity (alog the axis) [m/s2]
                                % horizontal: gx = 0
                                % vertical, going up: gx = -9.81 m/s2
                                % vertical, going down: gx =  9.81 m/s2
                                
                                
% ------------------------------------------------------------------------%
% 3. Simulation
% ------------------------------------------------------------------------%
clear screen;
disp( sprintf('1. Isothermal, without pressure drop'));
disp( sprintf('2. Isothermal, with pressure drop'));
disp( sprintf('3. Non-Isothermal, without pressure drop'));
disp( sprintf(' '));
prompt = '   What type of simulation? ';
simulation = input(prompt);


if (simulation == 1)
    
    % Initial conditions
    y0 = [omega0, T0, P0];
    
    % ODE Solver
    [x,y] = ode15s(@odeTubularIsothermal,[x0 xF], y0);
    
elseif (simulation == 2)

    % Initial conditions
    y0 = [omega0, T0, P0];
    
    % ODE Solver
    [x,y] = ode15s(@odeTubularIsothermalWithPressureDrop,[x0 xF], y0);
    
elseif (simulation == 3)
    
    % Initial conditions
    y0 = [omega0, T0, P0];
    
    % ODE Solver
    [x,y] = ode15s(@odeTubularHeatExchange,[x0 xF], y0);
    
end

% ------------------------------------------------------------------------%
% 3. Post processing
% ------------------------------------------------------------------------%

figure; % create new figure
hold all;

% Mass fraction profiles
for i=1:kinetics.ns
    subplot(2,2,1);
    hold all;
    plot (x, y(:,i),'LineWidth',2);
end
title('Mass fractions of species along the reactor');
xlabel('axial length [m]'); 
ylabel('mass fractions [-]'); 
legend(kinetics.species);
axis([-inf,inf,0,inf]);

% Conversion profile
subplot(2,2,2);
plot (x, 1-y(:,1),'LineWidth',2);
title('Conversion profile along the reactor');
xlabel('axial length [m]'); 
ylabel('Conversion');

% Temperature profile
subplot(2,2,3);
plot (x, y(:,kinetics.ns+1),'LineWidth',2);
title('Temperature profile along the reactor');
xlabel('axial length [m]'); 
ylabel('Temperature [K]');

% Pressure profile
subplot(2,2,4);
plot (x, y(:,kinetics.ns+2)/1e5,'LineWidth',2);
title('Pressure profile along the reactor');
xlabel('axial length [m]'); 
ylabel('Pressure [bar]');
