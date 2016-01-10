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
%   Copyright(C) 2016 Alberto Cuoci                                       |
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

function funs = PackedBedReactorFunctions
  funs.withpressuredrop=@withpressuredrop;
  funs.withoutpressuredrop = @withoutpressuredrop;
end

%*************************************************************************%
% Part 1: pure plug-flow reactor
%*************************************************************************%

function dX_over_dW = withpressuredrop(W,X)
    
    global k;
    global Kc;
    global Ka;
    global PAin;
    global FAin;
    global alfa;
    
    % Calculate the pressure drop
    y = sqrt(1-alfa*W);

    % Partial pressures of species    
    PA   = PAin*(1-X)*y;
    PB   = PAin*(1.5-X)*y;
    PC   = PAin*X*y;
    
    % Overal reaction rate
    rprime = k*PA*PB/(1+Kc*PC+Ka*PA);
    
    %Equation
    dX_over_dW = rprime/FAin;
    
end

%*************************************************************************%
% Part 2: membrane reactor (only B component)
%*************************************************************************%

function dX_over_dW = withoutpressuredrop(W,X)
    
    global k;
    global Kc;
    global Ka;
    global PAin;
    global FAin;
    
    % No pressure drop
    y = 1;

    % Partial pressures of species
    PA   = PAin*(1-X)*y;
    PB   = PAin*(1.5-X)*y;
    PC   = PAin*X*y;
    
    % Overal reaction rate
    rprime = k*PA*PB/(1+Kc*PC+Ka*PA);
    
    %Equation
    dX_over_dW = rprime/FAin;
    
end