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
%   Copyright(C) 2021 Alberto Cuoci                                       |
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
%                                                                         %
% CSTR (constant density)                                                 %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all;  clear variables;


%% Input data
k1 = 0.5;       % [1/min]
k2 = 0.1;       % [1/min]
CAin = 2;       % [mol/l]
CBin = 0;       % [mol/l]
CCin = 0;       % [mol/l]


%% Non linear system solution
nls_opt = optimset('Display','off');
for i=1:60
    
    Tau = i;                                    % [min]
    CFirstGuess = [CAin/3 CAin/3 CAin/3]';      % [mol/l]
    C = fsolve(@isothermalCSTR, CFirstGuess, nls_opt, k1,k2,CAin,CBin,CCin,Tau);
    
    t(i) = Tau;                     % [min]
    X(i) = (CAin-C(1))/CAin;        % [-]
    sigmaB(i) = C(2)/(CAin-C(1));   % [-]
    yieldB(i) = C(2)/CAin;          % [-]
    
end


%% Concentrations
figure;
plot(t,sigmaB,'-o', t,yieldB,'-o', t,X,'-o');
legend('sigmaB', 'yieldB', 'X');
xlabel('time [min]');
ylabel('sigma, yield, X');
title('Main results');


%% CSTR NLS (defined locally, requires at least MATLAB-2016b)
function F = isothermalCSTR(C, k1,k2,CAin,CBin,CCin,Tau)

    % Recover main variables
    CA = C(1);   % [mol/l]
    CB = C(2);   % [mol/l]
    CC = C(3);   % [mol/l]
    
    % Formation rates
    RA = -(k1+k2)*CA;   % [mol/l/min]
    RB = k1*CA;         % [mol/l/min]
    RC = k2*CA;         % [mol/l/min]

    % Mass balance equations
    F(1) = CAin-CA + RA*Tau;
    F(2) = CBin-CB + RB*Tau;
    F(3) = CCin-CC + RC*Tau;
    
end
