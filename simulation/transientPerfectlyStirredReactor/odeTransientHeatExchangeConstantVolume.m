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

% Isothermal transient perfectly stirred reactor
% Assumptions:
% 1. constant volume (rigid vessel)
% 2. constant temperature (isothermal conditions)
% 3. single inlet stream
% 4. no mass accumulation (i.e. min = mout)

function dy = odeTransientHeatExchangeConstantVolume(t,y)

global kinetics;
global omega_in;
global T_in;
global V0;
global mfr_in;
global Te;
global D;
global U;
global Le;

% 1. Pre-allocate vectors
dy = zeros(kinetics.ns+2,1);
omega = zeros(kinetics.ns,1);
domega =zeros(kinetics.ns,1);

% 2. Recover the variables
for i=1:kinetics.ns
    omega(i) = y(i);
end
T = y(kinetics.ns+1);
m = y(kinetics.ns+2);

% 3. Residence time [s]
Tau = m/mfr_in;

% 4. Density (definition) [kg/m3]
rho = m/V0;

% 5. Pressure [Pa]
P = kinetics.Pressure(rho,T,omega);

% 6. Reaction rates and formation rates [kmol/m3/s]
[r,R] = kinetics.Calculate(T,P,omega);

% ------------------------------------------------------------------------%
% 7. Species equations (differential)  [1/s]
% ------------------------------------------------------------------------%
for i=1:kinetics.ns
    domega(i) = (omega_in(i)-omega(i))/Tau +R(i)*kinetics.MW(i)/rho;
end

% ------------------------------------------------------------------------%
% 8. Energy equation (differential)  [K/m]
% ------------------------------------------------------------------------%
Cp = kinetics.SpecificHeat(omega);  % [J/kg/K]
Qr = kinetics.ReactionHeat(T,r,R);  % [W/m3]
h_in = kinetics.Enthalpies(T_in);   % [J/kg]
h_out= kinetics.Enthalpies(T);      % [J/kg]

delta_h = 0;
for i=1:kinetics.ns
    delta_h = delta_h + omega_in(i)*(h_in(i)-h_out(i));
end

% Sum of formation rate of species [kmol/m3/s]
sumR = 0.;
for i=1:kinetics.ns
   sumR = sumR+R(i);
end

% Beta coefficient and molecular weight
Beta = kinetics.Beta(T,P,omega);
W = kinetics.MolecularWeight(omega);

dT =    (Beta*T/rho*P/rho)*W*sumR + ...
        delta_h/Tau + Qr/rho + ...
		6/D*U*(Te-T)/rho + ...
		Le/(rho*V0);

dT = dT / (Cp - Beta*P/rho);

% ------------------------------------------------------------------------%
% 8. Recover equations
% ------------------------------------------------------------------------%
for i=1:kinetics.ns
    dy(i) = domega(i);
end
dy(kinetics.ns+1) = dT;
dy(kinetics.ns+2) = 0;
