%-------------------------------------------------------------------------%
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
% Example of maximum mixedness model                                      %
% Experimental RTD function available                                     %
% Non linear fitting of RTD                                               %
% Numerical solution of a ODE problem with IC                             %
% Requires at least MATLAB-2016b                                          %
%                                                                         %
%-------------------------------------------------------------------------%
close all; clear variables;

global k;            % [1/min]
global CAin;         % [kmol/m3]
global C3;

%% Experimental data
nExp = 10;
tExp = [0 0.5 1 1.5 2 2.5 3 3.5 4 5]';
FExp = [0.0421434 0.4344916 0.6686519 0.8122302 0.8811733 ...
        0.9367312 0.9622991 0.9786636 0.9827176 0.995]';

% User data
V    = 0.7;    % total volume [l]   
Q    = 0.7;    % total volumetric flow rate [l/min]
CAin = 10;     % [kmol/m3]
k    = 5;      % kinetic constant [m3/kmol/min]
tau  = V/Q;    % nominal residence time [min]

% Ideal CSTR and PFR solutions
CAcstr = (-1+sqrt(1+4*k*tau*CAin))/(2*k*tau);
CApfr = CAin/(1+k*tau*CAin);


%% Fitting of F using non-linear regression analysis
% This was already done in Ex. 3 using linear regression analysis
% Here we simply try an alternative method, but results are the same!

firstGuess = [1 tau];
nlm = fitnlm(tExp,FExp,@nonLinearModel,firstGuess);
C2   = nlm.Coefficients.Estimate(1);
C3   = nlm.Coefficients.Estimate(2);

% Quality evaluation
FModel = 1-C2*exp(-tExp/C3);
SSres = sum( (FExp-FModel).^2  );
SStot = sum( (FExp-mean(FExp)).^2 );
R2 = 1 - SSres/SStot;

% Residence Time Distribution function
EExp = C2/C3*exp(-tExp/C3);

% Print on the screen 
fprintf('C2: %f \n', C2);
fprintf('C3: %f \n', C3);
fprintf('R2 coefficient: %f \n', R2);


%% ODE solution
lambdaMax = 100*tau;
[lambda, CA] = ode23s(@ODEEquation, [lambdaMax 0], CAin);
CAout = CA(end);


%% Comparison
fprintf('A) Max. Mixedness: %f - CSTR: %f - PFR: %f [kmol/m3]\n', ...
         CAout, CAcstr, CApfr);
fprintf('X) Max. Mixedness: %f - CSTR: %f - PFR: %f \n', ...
         1-CAout/CAin, 1-CAcstr/CAin, 1-CApfr/CAin);
     

%% ODE system corresponding to the MM model
function dCA = ODEEquation(lambda,CA)

    global C3 k CAin;

    E_over_1_minus_F = 1/C3;
     
    RA = -k*CA^2;
    dCA = E_over_1_minus_F*(CA-CAin)-RA;

end


%% Function fitting the experimental CDF
function F = nonLinearModel(a,t)

    F = 1 - a(1)*exp(-t/a(2));

end