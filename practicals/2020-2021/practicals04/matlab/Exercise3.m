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
close all; clear variables;

%% Experimental data
Nexp = 9;
T =     [ 300        350        350        400        450        450        500        500        550        ]; % [K]
kappa = [ 2.5197e-04 3.4984e-02 3.8545e-02 1.5634e+00 2.9993e+01 2.5964e+01 3.3637e+02 2.9365e+02 2.1903e+03 ]; % [1/min]
      

%% Hypothesis: kappa=A*T^n*exp(-E/RT)
% Linearization: y = a0+a1*x1+a2*x2
%                y = ln(k), a0=ln(A), a1=n, a2=-E/R, x1=ln(T), x2=1/T

% Linear model
Yexp = log(kappa)';
X = [ ones(Nexp,1) log(T)' 1./T'];

% Linear system solution
% (X'*X)a = X'*Yexp
a = (X'*X)\(X'*Yexp);

% Recover original parameters
A = exp(a(1));
n = a(2);
E = -1.987*a(3);

% Statistics
Ymod = X*a;
SSres = (Yexp-Ymod)'*(Yexp-Ymod);
SStot = (Yexp-mean(Yexp))'*(Yexp-mean(Yexp));
R2 = 1 - SSres/SStot;

% Plot
figure; hold on;
plot( 1./T, log(kappa), '.','markersize',20);
plot( 1./T, log(A*T.^n.*exp(-E/1.987./T)), '-');
xlabel('1/T (1/K)'); ylabel('ln(k) (1/min)')

% Print on the screen 
fprintf('Frequency factor (1/min):    %e \n', A);
fprintf('Temperature exponent:        %f \n', n);
fprintf('Activation energy (cal/mol): %f \n', E);
fprintf('R2 coefficient:              %f \n', R2);
