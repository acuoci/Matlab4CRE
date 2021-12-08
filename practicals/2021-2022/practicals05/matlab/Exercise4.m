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
close all; clear variables;

%% Experimental data
Nexp = 7;
CA0 = [ 0.1;      0.3;     0.5;     0.8;     1;       2;       4;        ];
rA  = [ 0.073e-2; 0.32e-2; 0.77e-2; 1.43e-2; 1.8e-2;  4.7e-2;  12.34e-2; ];


%% Plot of raw data
figure; title('Experimental data');
scatter(CA0,rA, 'o'); 
xlabel('C_{A0} [mol/l]'); ylabel('r_A [mol/l/min]');


%% Hypothesis: rA = k*CA^n
% Linearization: y = q + m*x
%                y = ln(r), q = ln(k), m = n, x = ln(CA)
y = log(rA);    
x = log(CA0);
figure; title('Linearized model');
scatter(x,y, 'o'); 
xlabel('ln(C_{A0})'); ylabel('ln(r_A)');


%% Linear regression analysis
% Yexp = a0 + a1*X
% Yexp = ln(r), a0 = ln(k), a1 = n, X = [1, ln(CA0)]
Yexp = log(rA);
lnCA0 = log(CA0);
X = [ ones(Nexp,1) lnCA0 ];

% Linear system solution
% (X'*X)a = X'*Yexp
a = (X'*X)\(X'*Yexp);

% Recover original parameters
k = exp(a(1));
n = a(2);

% Statistics
Ymod = X*a;
SSres = (Yexp-Ymod)'*(Yexp-Ymod);
SStot = (Yexp-mean(Yexp))'*(Yexp-mean(Yexp));
R2 = 1 - SSres/SStot;

% Plot
figure;
plot( lnCA0, Yexp', '.','markersize',20);
refline(n, log(k));
xlabel('ln(C_{A0})');
ylabel('ln(r)');

% Print on the screen 
fprintf('Kinetic constant: %e \n', k);
fprintf('Reaction order:   %f \n', n);
fprintf('R2 coefficient:   %f \n', R2);
