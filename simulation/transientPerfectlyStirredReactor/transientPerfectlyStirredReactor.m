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
global omega_in;            % inlet mass fractions [-]
global mfr_in;              % mass flow rate (inlet) [kg/s]
global T_in;                % inlet temperature [K]
global V0;                  % initial volume [m3]
global P0;                  % initial pressure [Pa]
global D;                   % reactor diameter (pherical shape is assumed)
global U;                   % global heat exchange coefficient [W/m2/K]
global Te;                  % external temperature [K]
global Le;                  % power supplied by the mixer [W];

% ------------------------------------------------------------------------%
% 1. Load the kinetic mechanism 
% ------------------------------------------------------------------------%
kinetics = GasMechanism_AceticAnhydride;

% ------------------------------------------------------------------------%
% 2.a Define inlet stream
% ------------------------------------------------------------------------%
T_in = 1000.;                      % inlet temperature [K]
mfr_in = 0.1                       % inlet mass flow rate [kg/s]

omega_in = zeros(1,kinetics.ns);   % inlet mass fractions [-]
omega_in(1) = 1.;

% ------------------------------------------------------------------------%
% 2.b Define initial conditions
% ------------------------------------------------------------------------%
V0     = 1;                     % initial volume [m3]
P0     = 101325.;               % initial pressure [Pa]
T0     = 1000.;                 % initial temperature [K]
omega0 = zeros(1,kinetics.ns);  % initial mass fractions [-]
omega0(1) = 1.;

% ------------------------------------------------------------------------%
% 2.c Define additional data
% ------------------------------------------------------------------------%
D = (6*V0/pi)^(1/3);       % reactor diameter (pherical shape is assumed)
U = 200;                   % global heat exchange coefficient [W/m2/K]
Te = 1200;                 % external temperature [K]
Le = 1000;                 % power supplied by the mixer [W]

% ------------------------------------------------------------------------%
% 2.d Define integration range
% ------------------------------------------------------------------------%
t0 = 0;                         % initial time [s]
tF = 20;                        % final time [s]
                                
% ------------------------------------------------------------------------%
% 4. Simulation
% ------------------------------------------------------------------------%
clear screen;
disp( sprintf('1. Isothermal, constant volume'));
disp( sprintf('2. Isothermal, constant pressure'));
disp( sprintf('3. Non-Isothermal, constant volume'));
disp( sprintf(' '));
prompt = '   What type of simulation? ';
simulation = input(prompt);

if (simulation == 1)
    
    % Initial conditions
    y0 = [omega0, T0, kinetics.Density(T0,P0,omega0)*V0];
    
    % ODE Solver
    [t,y] = ode15s(@odeTransientIsothermalConstantVolume,[t0 tF], y0);
    
elseif (simulation == 2)
    
    % Initial conditions
    y0 = [omega0, T0, kinetics.Density(T0,P0,omega0)*V0];
    
    % ODE Solver
    [t,y] = ode15s(@odeTransientIsothermalConstantPressure,[t0 tF], y0);
    
elseif (simulation == 3)
        
    % Initial conditions
    y0 = [omega0, T0, kinetics.Density(T0,P0,omega0)*V0];
    
    % ODE Solver
    [t,y] = ode15s(@odeTransientHeatExchangeConstantVolume,[t0 tF], y0);
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
    plot (t, y(:,i),'LineWidth',2);
end
title('Mass fractions of species');
xlabel('time [s]'); 
ylabel('mass fractions [-]'); 
legend(kinetics.species);
axis([-inf,inf,0,inf]);

% Conversion profile
subplot(2,2,2);
plot (t, 1-y(:,1),'LineWidth',2);
title('Conversion profile');
xlabel('time [s]');
ylabel('Conversion');

% Temperature profile
subplot(2,2,3);
plot (t, y(:,kinetics.ns+1),'LineWidth',2);
title('Temperature profile');
xlabel('time [s]');
ylabel('Temperature [K]');

if (simulation == 1 || simulation == 3)
    % Pressure (post-processing of solution is needed)
    omega = zeros(1,kinetics.ns); 
    for i=1:size(t)
        for j=1:kinetics.ns
            omega(j) = y(i,j);
        end
        P(i) = kinetics.Pressure(y(i,kinetics.ns+2)/V0,y(i,kinetics.ns+1),omega);
    end

    subplot(2,2,4);
    plot (t, P/1e5,'LineWidth',2);
    title('Pressure');
    xlabel('time [s]');
    ylabel('pressure [bar]');
end

if (simulation == 2)
    % Volume (post-processing of solution is needed)
    omega = zeros(1,kinetics.ns); 
    for i=1:size(t)
        for j=1:kinetics.ns
            omega(j) = y(i,j);
        end
        V(i) = y(i,kinetics.ns+2)/kinetics.Density(y(i,kinetics.ns+1), P0,omega);
    end

    subplot(2,2,4);
    plot (t, V,'LineWidth',2);
    title('Volume');
    xlabel('time [s]');
    ylabel('volume [m3]');
end