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

function funs = MembraneReactorFunctions
  funs.plugflow=@plugflow;
  funs.membranereactor_only_B = @membranereactor_only_B;
end

%*************************************************************************%
% Part 1: pure plug-flow reactor
%*************************************************************************%

function dy = plugflow(t,y)
    
    global k;
    global KC;
    global T0;
    global P0;
    
    FA = y(1);                % molar flow rate of species A [mol/min]
    FB = y(2);                % molar flow rate of species A [mol/min]
    FC = y(3);                % molar flow rate of species A [mol/min]
    
    Ftot = FA+FB+FC;
    
    Ctot = P0/8314/T0;
    CA = Ctot*FA/Ftot;
    CB = Ctot*FB/Ftot;
    CC = Ctot*FC/Ftot;
    
    R  = k*(CA-CB*CC/KC);             % reaction rate [mol/dm3/min]
    
    dy = zeros(3,1);
    dy(1) = -R;                       % dFA_over_dV [mol/dm3/min]
    dy(2) =  R;                       % dFB_over_dV [mol/dm3/min]
    dy(3) =  R;                       % dFC_over_dV [mol/dm3/min]
    
end

%*************************************************************************%
% Part 2: membrane reactor (only B component)
%*************************************************************************%

function dy = membranereactor_only_B(t,y)
    
    global k;
    global KC;
    global T0;
    global P0;
    global km;    
    
    FA = y(1);                % molar flow rate of species A [mol/min]
    FB = y(2);                % molar flow rate of species A [mol/min]
    FC = y(3);                % molar flow rate of species A [mol/min]
    
    Ftot = FA+FB+FC;
    
    Ctot = P0/8314/T0;
    CA = Ctot*FA/Ftot;
    CB = Ctot*FB/Ftot;
    CC = Ctot*FC/Ftot;
    
    R  = k*(CA-CB*CC/KC);             % reaction rate [mol/dm3/min]
    
    dy = zeros(3,1);
    dy(1) = -R;                       % dFA_over_dV [mol/dm3/min]
    dy(2) =  R - km*CB;               % dFB_over_dV [mol/dm3/min]
    dy(3) =  R;                       % dFC_over_dV [mol/dm3/min]
    
end
