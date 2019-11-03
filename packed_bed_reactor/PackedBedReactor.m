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
% Packed Bed Reactor with pressure drop:                                  %
% For details see: Fogler, Chapter 10, Example 10.3                       %
%                                                                         %
%-------------------------------------------------------------------------%
close all; clear variables;


%% Inputa data

% Operating conditions
alpha = 9.8e-5;   % pressure drop coefficient [1/kg]
FAin  = 50;       % inlet molar flow rate of A [mol/min]
Pin   = 40;       % pressure [atm]
T     = 640;      % temperature [C]
rhoB  = 400;      % bed density [kg/m3]
xAin  = 0.30;     % inlet molar fraction of A
xBin  = 0.45;     % inlet molar fraction of B
xIin  = 0.25;     % inlet molar fraction of inerts

% Kinetic and thermodynamic data
k    = 0.00087;   % kinetic parameter [mol/atm2/kgcat/min]
Ka   = 1.39;      % kinetic parameter [1/atm]
Kc   = 1.038;     % kinetic parameter [1/atm] 

% Maximum amount of catalyst
Pmin = 1;                        % minimum (relative pressure) [atm]
Wmax = (1-(Pmin/Pin)^2) / alpha; % maximum amount of catalyst [kg]

% Initial conditions
PAin = Pin*xAin;    % inlet partial pressure of A [atm]
X0   = 0.;          % initial conversion of species A


%% Part 1: PBR with pressure drop

% ODE solution
[W,X] = ode45(@withpressuredrop, [0 0.98*Wmax], X0, [], k,Kc,Ka,PAin,FAin,alpha);

% Plotting solution
y  = sqrt(1-alpha*W);
PA = PAin*(1-X).*y;
PB = PAin*(1.5-X).*y;
PC = PAin*X.*y;
    
figure;
subplot(2,3,1);
hold all;
plot(W, PA, 'LineWidth', 2);
plot(W, PB, 'LineWidth', 2);
plot(W, PC, 'LineWidth', 2, 'LineStyle', '--');
title ('PBR with pressure drop');   legend('A','B','C');
xlabel('mass of catalyst [kg]');    ylabel('pressures [atm]');

subplot(2,3,2);
hold all;
plot(W,Pin*y, 'LineWidth', 2);
title ('PBR with pressure drop');   legend('Pressure');
xlabel('mass of catalyst [kg]');    ylabel('total pressure [atm]');

subplot(2,3,3);
hold all;
plot(W,X, 'LineWidth', 2);
title ('PBR with pressure drop');   legend('Conversion');
xlabel('mass of catalyst [kg]');    ylabel('conversion');


%% Part 2: PBR without pressure drop

% ODE solution
[W,X] = ode45(@withoutpressuredrop, [0 0.98*Wmax], X0, [], k,Kc,Ka,PAin,FAin);

% Plotting solution
PA = PAin*(1-X);
PB = PAin*(1.5-X);
PC = PAin*X;

subplot(2,3,4);
hold all;
plot(W, PA, 'LineWidth', 2);
plot(W, PB, 'LineWidth', 2);
plot(W, PC, 'LineWidth', 2, 'LineStyle', '--');
title ('PBR without pressure drop');    legend('A','B','C');
xlabel('mass of catalyst [kg]');    ylabel('pressures [atm]');

subplot(2,3,5);
hold all;
plot(W,Pin.*W./W, 'LineWidth', 2);
title ('PBR without pressure drop');    legend('Pressure');
xlabel('mass of catalyst [kg]');    ylabel('total pressure [atm]');
axis([-Inf Inf Pin*0.99 Pin*1.01]);

subplot(2,3,6);
hold all;
plot(W,X, 'LineWidth', 2);
title ('PBR without pressure drop');    legend('Conversion');
xlabel('mass of catalyst [kg]');    ylabel('conversion');



%% ODE system 1: PBR with pressure drop
function dX_over_dW = withpressuredrop(W,X, k,Kc,Ka,PAin,FAin,alpha)
    
    % Calculate the pressure drop
    y = sqrt(1-alpha*W);

    % Partial pressures of species    
    PA   = PAin*(1-X)*y;
    PB   = PAin*(1.5-X)*y;
    PC   = PAin*X*y;
    
    % Overal reaction rate
    rprime = k*PA*PB/(1+Kc*PC+Ka*PA);
    
    %Equation
    dX_over_dW = rprime/FAin;
    
end


%% ODE system 2: PBR without pressure drop
function dX_over_dW = withoutpressuredrop(W,X, k,Kc,Ka,PAin,FAin)
    
    % No pressure drop
    y = 1;

    % Partial pressures of species
    PA   = PAin*(1-X)*y;
    PB   = PAin*(1.5-X)*y;
    PC   = PAin*X*y;
    
    % Overal reaction rate
    rprime = k*PA*PB/(1+Kc*PC+Ka*PA);
    
    %Equation
    dX_over_dW = rprime/FAin;
    
end
