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

