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

% The residence times for a single, user-defined reaction
% in different reactors: batch (BR), plug flow (FPR), and CSTR:
%
% The residence times are calculated through numerical integration of
% design equations. The reaction rates as a function of the conversion are
% specified by the user through the RA(X) function

close all;
clear all;

% ------------------------------------------------------------------------%
% Initial concentration [kmol/m3/s]
% ------------------------------------------------------------------------%
CA0 = 1;

% ------------------------------------------------------------------------%
% User defined reaction rate function
% The reaction rates has to be returned in [kmol/me/s]
% ------------------------------------------------------------------------%

% Example: irreversible, second order reaction
% ------------------------------------------------------------------------%
% A + B -> 2C
% r = k1*CA*CB
% ------------------------------------------------------------------------%
k1  = 1;        % kinetic constant [m3/kmol/s]
eps = 0;        % changes of moles
tetaB0 = 2;     % relative amount of initial B
RA = @(x) -k1* ( CA0*(1-x)./(1+eps*x) ) .* ( CA0*(tetaB0-x)./(1+eps*x) )  ;

% Example: equilibrium reaction
% ------------------------------------------------------------------------%
% A <-> B
% r = k1*CA+k1/KC*CB
% ------------------------------------------------------------------------%
%k1  = 1;       % kinetic constant [1/s]
%KC = 10;       % equilibrium constant
%tetaB0 = 0;    % relative amount of initial B
%RA = @(x) -k1*CA0./(1+eps*x).*(1-x-(tetaB0+x)/KC);


% ------------------------------------------------------------------------%
% Design equations for ideal reactors (functions to be integrated)
% ------------------------------------------------------------------------%
integral_PFR  = @(x) CA0./(-RA(x));
integral_CSTR = @(x) CA0.*x/(-RA(x));


% ------------------------------------------------------------------------%
% Create plots
% ------------------------------------------------------------------------%
figure;     % create new figure
np=500;     % number of points (only for post-processing)
X = linspace(0,0.999,np);   % create the x axis for plotting profiles

% ------------------------------------------------------------------------%
% Integrate design equations
% ------------------------------------------------------------------------%
for i=1:np
TauPFR(i) = integral(@(x)integral_PFR(x),0,X(i));
TauCSTR(i) = integral_CSTR(X(i));
end

% ------------------------------------------------------------------------%
% Clean data for best graphical results
% ------------------------------------------------------------------------%
np_extreme = np;
for i=np:-1:2
    if (TauPFR(i)<=0)
        if (i<np_extreme)   
            np_extreme = i;
        end
        TauPFR(i) = 1e-20;
    end
    if (TauCSTR(i)<=0)
        if (i<np_extreme)   
            np_extreme = i;
        end
        TauCSTR(i) = 1e-20;
    end
end

% ------------------------------------------------------------------------%
% Plot
% ------------------------------------------------------------------------%
subplot(1,2,1);
loglog( 1-X(1:np_extreme),TauPFR(1:np_extreme), '.', ...
        1-X(1:np_extreme),TauCSTR(1:np_extreme), '.');
title('PFR vs CSTR');
%xlim([0.001 1]);
ylim([0.1 1000]);
xlabel('1-X');
ylabel('residence time [ut]');
legend('Plug Flow', 'CSTR');

%Plot
subplot(1,2,2);
loglog(1-X(1:np_extreme),TauCSTR(1:np_extreme)./TauPFR(1:np_extreme), '.');
title('PFR vs CSTR');
%xlim([0.001 1]);
ylim([0.1 1000]);
xlabel('1-X');
ylabel('TauCSTR/TauPFR');
