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

%Kinetic mechanism defining the following reactions

% CH3COCH3 => CH2CO + CH4 (r1 = k1*CA)

classdef GasMechanism_AceticAnhydride < KineticMechanism

    properties
        
    end
    
    methods
        
        function obj = GasMechanism_AceticAnhydride
        
        % number of reactions
        obj.nr = 1;
    
        % number of species
        obj.ns = 3;  
        
        % names of species
        obj.species = char('CH3COCH3','CH2CO','CH4');
    
        % molecular weights [kg/kmol]
        obj.MW = [58., 42., 16.];
    
        % formation enthalpy [kJ/mol]
        obj.H0 = [-216.67, -61.09, -74.81];
    
        % specific heats [J/mol/K]
        obj.Cp = [163., 83., 71.];
    
        % reference temperature [K]
        obj.Tref = 298;
        
        end
        
        % density [kg/m3]
        function rho = Density(obj, T,P,omega)
          
            rho = P*obj.MolecularWeight(omega)/obj.Rgas/T;
        
        end
        
        % pressure [Pa]
        function P = Pressure(obj, rho, T, omega)
          
            P = rho*obj.Rgas*T/obj.MolecularWeight(omega);
        
        end
        
        % beta coefficient [1/T]
        function beta = Beta(obj, T,P,omega)
          
            beta = 1/T;
        
        end
        
        % viscosity [kg/m/s]
        function mu = Viscosity(obj, T,P,omega)
          
           mu = 1.85e-5;
        
        end
        
        % reaction rates [kmol/m3/s]
        function r = ReactionRates(obj, T,P,C)
          
            r = zeros(obj.nr,1);
            r(1) = 3.58*exp(34222*(1/1035-1/T))*C(1);
        
        end
        
        % formation rates [kmol/m3/s]
        function R = FormationRates(obj, r)
          
            R = zeros(obj.ns,1);
            
            R(1) = -r(1);
            R(2) =  r(1);
            R(3) =  r(1);
        
        end
        
        
    end
end


