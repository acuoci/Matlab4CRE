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
% Example of compartment model                                            %
% Bypass + dead volume (2 parameter model)                                %
% Parameters from linear regressions analysis                             %
% Requires at least MATLAB-2016b                                          %
%-------------------------------------------------------------------------%
close all; clear variables;

global k;            % [1/min]
global CAin;         % [kmol/m3]
global tauEff;       % [min]
global alpha;        % [-]


%% Experimental data
nExp = 10;
tExp = [0 0.5 1 1.5 2 2.5 3 3.5 4 5]';
FExp = [0.0421434 0.4344916 0.6686519 0.8122302 0.8811733 ...
        0.9367312 0.9622991 0.9786636 0.9827176 0.995]';

% User data
V    = 0.7;    % total volume [l]   
Q    = 0.7;    % total volumetric flow rate [l/min]
CAin = 10;     % inlet concentration of A [kmol/m3]
k    = 5;      % kinetic constant [m3/kmol/min]
tau  = V/Q;    % nominal residence time [min]

% Ideal CSTR solution
CAcstr = (-1+sqrt(1+4*k*tau*CAin))/(2*k*tau);


%% Analysis of experimental data

% Linear regression
Y = log(1./(1-FExp));
X = [ones(nExp,1) tExp];
a = (X'*X)\(X'*Y);

% Statistics
SSres = (Y-X*a)'*(Y-X*a);
SStot = (Y-mean(Y))'*(Y-mean(Y));
R2 = 1 - SSres/SStot;

% Calculation of compartment model parameters
q = a(1);
m = a(2);
alpha = 1-1/exp(q);
beta  = 1-(1-alpha)/m/tau;
Veff = (1-beta)*V;
Qeff = (1-alpha)*Q;
tauEff = Veff/Qeff;

% Print on the screen 
fprintf('alpha: %f \n', alpha);
fprintf('beta:  %f \n', beta);
fprintf('R2 coefficient: %f \n', R2);


%% NLS solution
YFirstGuess = [CAin CAin]';
options = optimset('Display', 'off');
Y = fsolve(@CSTR, YFirstGuess, options);
CA1 = Y(1);  
CAout = Y(2);

% Comparison
fprintf('A) Compartment model: %f - CSTR: %f [kmol/m3]\n', CAout, CAcstr);
fprintf('X) Compartment model: %f - CSTR: %f [kmol/m3]\n', 1-CAout/CAin, 1-CAcstr/CAin);


%% Non linear System (NLS)
function F = CSTR(C)

    global k;            % [m3/kmol/min]
    global CAin;         % [kmol/m3]
    global tauEff;       % [min]
    global alpha;        % [-]
    
    % Recover variables
    CA1   = C(1);
    CAout = C(2);
    
    % Reaction rates
    r = k*CA1^2;           % [kmol/m3/min]
    
    % Formation rates
    RA = -r;              % [kmol/m3/min]
    
    % Non linear equations
    F(1) = (CAin-CA1)/tauEff + RA;
    F(2) = CA1*(1-alpha) + CAin*alpha - CAout;
    
end