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

% Isothermal transient perfectly stirred reactor
% Assumptions:
% 1. constant volume (rigid vessel)
% 2. constant temperature (isothermal conditions)
% 3. single inlet stream
% 4. no mass accumulation (i.e. min = mout)

function dy = odeTransientIsothermalConstantVolume(t,y, kinetics,omega_in,V0,mfr_in)

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

    % 5. Pressure (from equation of state) [Pa]
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
    % 8. Recover equations
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        dy(i) = domega(i);
    end
    dy(kinetics.ns+1) = 0;
    dy(kinetics.ns+2) = 0;
    
end

