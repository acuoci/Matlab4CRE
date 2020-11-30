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
global CA0;

%% Experimental data
Nexp = 9;
CA =  [ 10 8.0262 6.5575 5.4393 4.5708 3.8847 3.3346 2.8877 2.5203 ];   % mol/l
t  =  [ 0  5      10     15     20     25     30     35     40     ];   % min


%% Approach 1: linear regression analysis (iterative procedure)
% Reaction orders to be tested
n = [1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60];

% Hypothesis: rA = k*CA^n
% Integration:   CA^(1-n) = CA0^(1-n) -k*(1-n)*t
% Linearization: y = q + m*x
%                y = CA^(1-n)-CA0^(1-n), q = 0, m = -k*(1-n), x = t

for j=1:length(n)
    
    % Yexp = a0 + a1*X
    % Yexp = CA^(1-n)-CA0^(1-n), a0 = 0, a1 = -k*(1-n), X = [t]
    CA0 = CA(1);
    Yexp = (CA.^(1-n(j))-CA0^(1-n(j)))';
    Xexp = [ t' ];

    % Linear system solution
    % (X'*X)a = X'*Yexp
    a = (Xexp'*Xexp)\(Xexp'*Yexp);

    % Recover original parameters
    k = -a(1)/(1-n(j));
    
    % Statistics
    Ymod = Xexp*a;
    SSres = (Yexp-Ymod)'*(Yexp-Ymod);
    SStot = (Yexp-mean(Yexp))'*(Yexp-mean(Yexp));
    R2 = 1 - SSres/SStot;
    
    fprintf('n=%f, k=%f, R2=%f\n', n(j), k, R2);
    
end


%% Approach 2: non-linear regression analysis

% Non linear regression
firstGuess = [1.e-2 1.5];
X = [CA' t'];
Y = zeros(Nexp,1);
nlm = fitnlm(X,Y,@nonLinearModel,firstGuess);

% Results
k = nlm.Coefficients.Estimate(1);
n = nlm.Coefficients.Estimate(2);
fprintf('Non linear regressions analaysis\n');
fprintf('n=%f, k=%f\n', n, k);

%% Function to be regressed
function y = nonLinearModel(a,x)

    global CA0;
    
    k     = a(1);
    n     = a(2);
   
    CA = x(:,1);
    t  = x(:,2);
    
    y = CA.^(-n+1) - CA0^(-n+1) +k*(1-n)*t;

end
