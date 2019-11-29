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
close all; clear variables;

%% Experimental data
nExp = 18;
PCO  = [ 1 1.8 4.08 1   1   1 2   2   2 1 1.8 4.08 3   3   3 1 1.8 4.08 ]';
PH2  = [ 1 1   1    0.1 0.5 4 0.1 0.5 4 2 2   2    0.1 0.5 4 3 3   3    ]';
rCH4 = [ 0.0072; 0.0129; 0.0292;  0.0049; 0.0073; 0.0053; 0.0098; ...
         0.0146; 0.0106; 0.0064; 0.0115;  0.0260; 0.0147; 0.0219; ...
         0.0159; 0.0058; 0.0104; 0.0235;];

     
%% Plot of raw data
figure;title('Experimental data: P_{CO}');
scatter(PCO, rCH4);
xlabel('P_{CO} [atm]'); ylabel('r [mol/kg/min]');

figure;title('Experimental data: P_{H2}');
scatter(PH2, rCH4);
xlabel('P_{H2} [atm]'); ylabel('r [mol/kg/min]');


%% Hypothesis
% r = b1*PCO^nCO*PH2^nH2low
%     ---------------------
%      1 + b2*PH2^nH2high

% Linearization (partial)
% y = log(b1) + nCO*log(PCO) + nH2low*log(PH2) - ...
%     log(1+b2*PH2.^nH2high);
% y = a0 + a1*x1 + a2*x2 + a3*x3
% a0 = log(b1), a1=nCO, a2=nH2low, a3=1
% x1 = log(PCO), x2 = log(PH2), x3 = - log(1+b2*PH2^nH2high)

Yexp = log(rCH4);
lnPCO0 = log(PCO);
lnPH20 = log(PH2);

eps_best = 1e10;
for i=1:10
    for j=1:10
        
        % Fix parameter values
        nH2high = 0.5 + 0.1*i;
        b2 = 1.5 + 0.1*j;
        
        % Linear regression
        X = [ ones(nExp,1) lnPCO0 lnPH20 -log(1+b2*PH2.^nH2high)];
        a = (X'*X)\(X'*Yexp);
        
        % Statistics (linear model)
        Ymod = X*a;
        SSres = (Yexp-Ymod)'*(Yexp-Ymod);
        SStot = (Yexp-mean(Yexp))'*(Yexp-mean(Yexp));
        R2 = 1 - SSres/SStot;
        
        % Check for the best model
        eps = abs(a(4)-1);
        if (eps<eps_best)
            eps_best = eps;
            best = [a(1),a(2),a(3),a(4), nH2high, b2, R2];
        end
        
    end
end

% Recover best parameters
a = [best(1) best(2) best(3) best(5) best(6)];
b1 = exp(a(1));
nCO = a(2);
nH2low = a(3);
nH2high = a(4);
b2 = a(5);
R2 = best(7);

% Print on the screen 
fprintf('Kinetic parameter nCO:     %e \n', nCO);
fprintf('Kinetic parameter b1:      %e \n', b1);
fprintf('Kinetic parameter b2:      %e \n', b2);
fprintf('Kinetic parameter nH2low:  %e \n', nH2low);
fprintf('Kinetic parameter nH2high: %e \n', nH2high);
fprintf('R2 coefficient:            %e \n', R2);
