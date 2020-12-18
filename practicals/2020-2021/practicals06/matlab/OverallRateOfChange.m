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
%   Copyright(C) 2020 Alberto Cuoci                                       |
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

function [rIIII,E,MH,resL,resG,resK] = OverallRateOfChange(pA, CB, data)

    %% Extract data from data struct
    b = data.b;             % stoichiometric coefficient species B
    GammaA = data.GammaA;   % diffusion coefficient [m2/s]
    GammaB = data.GammaB;   % diffusion coefficient [m2/s]
    fL = data.fL;           % liquid fraction [-]
    HA = data.HA;           % Henry's constant [Pa.m3/kmol]
    a = data.a;             % interfacial area per unti ov volume [m2/m3]
    KLa = data.KLa;         % mass transfer coefficient (liquid side) [1/s]
    KGa = data.KGa;         % mass transfer coefficient (gas side) [kmol/m3/s/Pa]
    kl = data.kl;           % kinetic constant (liquid volume basis) [m6/kmol^2/s]
    
    
    %% Preliminary calculations
    klStar  = kl*CB;                % [m3/kmol/s]
    klPrime = klStar*CB;            % [1/s]
    KL = KLa/a;                     % [m/s]
    MH=sqrt(klStar*GammaA*CB)/KL;   % [-]
    resG=1/(KGa);                   % [m3.s.Pa/kmol]
    resK=HA/(klPrime*fL);           % [m3.s.Pa/kmol]

    %% Enhancement factor
    if (MH<1)               % Regions I and II
        if (MH<0.3)
            E = 1;          % [-]   Region I
        else
            E = 1+MH^2/3;   % [-]   Region II
        end

        resL = HA/KLa/E;                % [m3.s.Pa/kmol]
        rIIII = pA/(resG+resL+resK);    % [kmol/m3/s]

    else % Regions III

        pAI = pA;                       % [Pa]
        for i=1:100
            Ei = 1+GammaB/GammaA*CB*HA/pAI/b;   % [-]

            if (MH<5*Ei)
                E = MH;
            else
                E = Ei;
            end

            resL = HA/KLa/E;                % [m3.s.Pa/kmol]
            rIIII = pA/(resG+resL+resK);    % [kmol/m3/s]

            rIIIIGas = KGa*(pA-pAI);        % [kmol/m3/s]
            if ( abs(rIIII-rIIIIGas)/rIIII<0.001)
                break;
            end

            pAI = pA-rIIII/(KGa);       % [Pa]
        end 
    end
    
end
