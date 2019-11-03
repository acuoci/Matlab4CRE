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

% Parent class defining a kinetic mechanism

classdef KineticMechanism

    properties
    
    % number of reactions
    nr;
    
    % number of species
    ns;
    
    %names of species
    species;
    
    % molecular weights [kg/kmol]
    MW;
    
    % formation enthalpy [kJ/mol]
    H0;
    
    % formation enthalpy [kJ/mol]
    dH0;    
    
    % specific heats [J/mol/K]
    Cp;
    
    % reference temperature [K]
    Tref;
    
    % ideal gas constant [J/kmol/K]
    Rgas = 8314.;
    
    end
    
    methods
        
        % molecular weight of the mixture
        function sum = MolecularWeight(obj, omega)
          
          sum = 0.;
          for i=1:obj.ns
              sum = sum + omega(i)/obj.MW(i);
          end
          
          sum = 1/sum;
          
        end
        
        % mole fractions from mass fraction
        function x = MoleFractions(obj, omega)
          
          MWmix = obj.MolecularWeight(omega);
          
          x = zeros(obj.ns,1);
          for i=1:obj.ns
              x(i)= omega(i)*MWmix/obj.MW(i);
          end
          
        end
        
        % calculates the reaction rates and the formation rates, both in [kmol/m3/s]
        function [r,R] = Calculate(obj, T, P, omega)
          
            x = obj.MoleFractions(omega);
            
            Ctot = obj.Density(T,P,omega)/obj.MolecularWeight(omega);
            C = Ctot*x;

            r = obj.ReactionRates(T,P,C);
            R = obj.FormationRates(r);
        
        end
        
        % calculates the mixture specific heat [J/kg/K]
        function CpMix = SpecificHeat(obj, omega)
          
          CpMix = 0;
          for i=1:obj.ns
              CpMix = CpMix + omega(i) * obj.Cp(i)*1.e3/obj.MW(i);
          end
          
        end
        
        % calculates the mixture enthalpy [J/kg]
        function HMix = Enthalpy(obj, T, omega)
          
          H = obj.Enthalpies(T);
          HMix = 0;
          for i=1:obj.ns
              HMix = Hmix + omega(i) * H(i)*1.e6/obj.MW(i);
          end
          
        end
        
        % calculates the enthalpies of the species [J/kg]
        function H = Enthalpies(obj, T)
          
          H = zeros(obj.ns,1);
          for i=1:obj.ns
              H(i) = obj.H0(i)*1.e6 + obj.Cp(i)*1.e3*(T-obj.Tref);
              H(i) = H(i)/obj.MW(i);
          end
          
        end
       
        % calculates the reaction heat [J/m3/s]
        function Qr = ReactionHeat(obj, T, r, R)

            if (size(obj.H0) ~= 0)
                
                H = obj.Enthalpies(T);    % [J/kg]
            
                Qr = 0;
                for i=1:obj.ns
                    Qr = Qr - R(i)*obj.MW(i)*H(i);
                end
                
            else

                Qr = 0;
                for i=1:obj.nr
                    Qr = Qr - r(i)*obj.dH0(i)*1e6;
                end
                
            end
            
        end
    end
end

