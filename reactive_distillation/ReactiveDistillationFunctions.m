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

function funs = ReactiveDistillationFunctions
  funs.semibatch=@semibatch;
  funs.reactive_distillation_only_C = @reactive_distillation_only_C;
  funs.reactive_distillation_every = @reactive_distillation_every;
end

%*************************************************************************%
% Part 1: pure semi-batch reactor (no reactive distillation)
%*************************************************************************%

function dy = semibatch(t,y)
    
    global k;
    global KC;
    global FBin;
    global Qin;

    V = y(5);                   % volume of liquid phase [dm3]
    CA = y(1)/V;                % concentration of species A [mol/dm3]
    CB = y(2)/V;                % concentration of species B [mol/dm3]
    CC = y(3)/V;                % concentration of species C [mol/dm3]
    CD = y(4)/V;                % concentration of species D [mol/dm3]
    
    R  = k*(CA*CB-CC*CD/KC);    % reaction rate [mol/dm3/min]
    
    dy = zeros(5,1);
    dy(1) = -R*V;               % dNA_over_dt [mol/min]
    dy(2) = -R*V + FBin;        % dNB_over_dt [mol/min]
    dy(3) = R*V;                % dNC_over_dt [mol/min]
    dy(4) = R*V;                % dND_over_dt [mol/min]
    dy(5) = Qin;                % dV_over_dt  [dm3/min]
    
end

%*************************************************************************%
% Part 2: semi-batch reactor + reactive distillation (only C component)
%*************************************************************************%

function dy = reactive_distillation_only_C(t,y)
    
    global k;
    global KC;
    global FBin;
    global Qin;
    global Ptot;
    global FAir; 
    global PvC;
    global CCout;

    V = y(5);                   % volume of liquid phase [dm3]
    CA = y(1)/V;                % concentration of species A [mol/dm3]
    CB = y(2)/V;                % concentration of species B [mol/dm3]
    CC = y(3)/V;                % concentration of species C [mol/dm3]
    CD = y(4)/V;                % concentration of species D [mol/dm3]
    
    R  = k*(CA*CB-CC*CD/KC);    % reaction rate [mol/dm3/min]
    
    XC = y(3)/(y(1)+y(2)+y(3)+y(4));  % mole fraction of C [-]
    FCout=XC*(PvC/Ptot)*FAir;         % evap. molar flow rate C [mol/min]
    
    dy = zeros(5,1);
    dy(1) = -R*V;                       % dNA_over_dt [mol/min]
    dy(2) = -R*V + FBin;                % dNB_over_dt [mol/min]
    dy(3) =  R*V - FCout;               % dNC_over_dt [mol/min]
    dy(4) =  R*V;                       % dND_over_dt [mol/min]
    dy(5) =  Qin - FCout/CCout;         % dV_over_dt  [dm3/min]
    
end

%*************************************************************************%
% Part 3: semi-batch reactor + reactive distillation (every component)
%*************************************************************************%

function dy = reactive_distillation_every(t,y)
    
    global k;
    global KC;
    global FBin;
    global Qin;
    global Ptot;
    global FAir; 

    global PvA;
    global PvB;
    global PvC;
    global PvD;

    global CAout;
    global CBout;
    global CCout;
    global CDout;

    V = y(5);                   % volume of liquid phase [dm3]
    CA = y(1)/V;                % concentration of species A [mol/dm3]
    CB = y(2)/V;                % concentration of species B [mol/dm3]
    CC = y(3)/V;                % concentration of species C [mol/dm3]
    CD = y(4)/V;                % concentration of species D [mol/dm3]
    
   R  = k*(CA*CB-CC*CD/KC);    % reaction rate [mol/dm3/min]
    
    Ntot = y(1)+y(2)+y(3)+y(4); % total number of moles [mol]
    XA = y(1)/Ntot;             % mole fraction of A [-]
    XB = y(2)/Ntot;             % mole fraction of B [-]
    XC = y(3)/Ntot;             % mole fraction of C [-]
    XD = y(4)/Ntot;             % mole fraction of D [-]
    
    FAout=XA*(PvA/Ptot)*FAir;   % evap. molar flow rate A [mol/min]
    FBout=XB*(PvB/Ptot)*FAir;   % evap. molar flow rate B [mol/min]
    FCout=XC*(PvC/Ptot)*FAir;   % evap. molar flow rate C [mol/min]
    FDout=XD*(PvD/Ptot)*FAir;   % evap. molar flow rate D [mol/min]
    
    dy = zeros(5,1);
    dy(1) = -R*V - FAout;               % dNA_over_dt [mol/min]
    dy(2) = -R*V - FBout + FBin;        % dNB_over_dt [mol/min]
    dy(3) =  R*V - FCout;               % dNC_over_dt [mol/min]
    dy(4) =  R*V - FDout;               % dND_over_dt [mol/min]
    dy(5) =  Qin - FAout/CAout ...
                 - FBout/CBout ...
                 - FCout/CCout ...
                 - FDout/CDout;         % dV_over_dt [dm3/min]
    
end