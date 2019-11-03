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

% Adiabatic tubular reactor
% Assumptions:
% 1. no pressure drop (low velocity and negligible viscous forces)
% 2. heat exchange (through a global exchange coefficient)
% 3. constant cross section area (i.e. constant internal diameter)

function dy = odeTubularHeatExchange(Tau,y, kinetics,T0,P0,v0,omega0,U,Te,D)

    % 1. Pre-allocate vectors
    dy = zeros(kinetics.ns+2,1);
    omega = zeros(kinetics.ns,1);
    domega =zeros(kinetics.ns,1);

    % 2. Recover the variables
    for i=1:kinetics.ns
        omega(i) = y(i);
    end
    T = y(kinetics.ns+1);
    P = y(kinetics.ns+2);

    % 3. Density (equation of state) [kg/m3]
    rho = kinetics.Density(T,P,omega);

    % 4. Velocity (continuity equation, algebraic) [m/s]
    rho0 = kinetics.Density(T0,P0,omega0);
    v   = rho0*v0/rho;

    % 5. Reaction rates and formation rates [kmol/m3/s]
    [r,R] = kinetics.Calculate(T,P,omega);

    % ------------------------------------------------------------------------%
    % 7. Species equations (differential)  [1/m]
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        domega(i) = R(i)*kinetics.MW(i)/(rho*v);
    end

    % ------------------------------------------------------------------------%
    % 8. Energy equation (differential)  [K/m]
    % ------------------------------------------------------------------------%
    Cp = kinetics.SpecificHeat(omega); % [J/kg/K]
    Qr = kinetics.ReactionHeat(T,r,R); % [W/m3]

    dT = Qr/(rho*Cp*v) + 4/D*U*(Te-T)/(rho*Cp*v);

    % ------------------------------------------------------------------------%
    % 9. Recover equations
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        dy(i) = domega(i);
    end
    dy(kinetics.ns+1) = dT;
    dy(kinetics.ns+2) = 0;

end
