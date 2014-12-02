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

% The residence times for a single, irreversible reaction are compared
% in different reactors: batch (BR), plug flow (FPR), and CSTR:
%
%       A -> Products      r = k*CA^n
%
% The residence times are calculated through numerical integration of
% design equations. The reaction rates as a function of the conversion are
% calculated using the analytical solution

close all;
clear all;

% ------------------------------------------------------------------------%
% Design equations for ideal reactors (functions to be integrated)
% ------------------------------------------------------------------------%
integral_BR = @(x,n,eps) (1+eps*x).^(n-1)./(1-x).^n;
integral_PFR = @(x,n,eps) (1+eps*x).^n./(1-x).^n;
integral_CSTR = @(x,n,eps) x*(1+eps*x).^n./(1-x).^n;

% ------------------------------------------------------------------------%
% Kinetic constant [kmol, m3, s] (the units depend on the reaction order n)
% ------------------------------------------------------------------------%
k1 = 1;

% ------------------------------------------------------------------------%
% Initial concentration [kmol/m3/s]
% ------------------------------------------------------------------------%
CA0 = 1;

% ------------------------------------------------------------------------%
% Order of reaction
% ------------------------------------------------------------------------%
prompt = 'What is the reaction order n? ';
n = input(prompt);

% ------------------------------------------------------------------------%
% Create plots
% ------------------------------------------------------------------------%
figure;                     % create new figure
X = linspace(1e-2,1-1e-2);  % create the x axis for plotting profiles

% ------------------------------------------------------------------------%
% 1. Case 1: eps = 0 (no changes in number of moles)
% ------------------------------------------------------------------------%
eps = 0;
for i=1:100
TauBR(i) = integral(@(x)integral_BR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauPFR(i) = integral(@(x)integral_PFR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauCSTR(i) = integral_CSTR(X(i),n,eps) / (k1*CA0^(n-1));
end

subplot(2,2,1);
loglog(1-X,TauBR, 1-X,TauPFR,'.', 1-X,TauCSTR);
title('eps = 0');
xlim([0.01 1]);
xlabel('1-X');
ylabel('residence time [s]');
legend('Batch','Plug Flow', 'CSTR');

% ------------------------------------------------------------------------%
% 2. Case 2: eps = 1
% ------------------------------------------------------------------------%
eps = 1;
for i=1:100
TauBR(i) = integral(@(x)integral_BR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauPFR(i) = integral(@(x)integral_PFR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauCSTR(i) = integral_CSTR(X(i),n,eps) / (k1*CA0^(n-1));
end

subplot(2,2,2);
loglog(1-X,TauBR, 1-X,TauPFR,'.', 1-X,TauCSTR);
title('eps = 1');
xlim([0.01 1]);
xlabel('1-X');
ylabel('residence time [s]');
legend('Batch','Plug Flow', 'CSTR');

% ------------------------------------------------------------------------%
% 3. Case 3: eps = -0.5
% ------------------------------------------------------------------------%
eps = -0.5;
for i=1:100
TauBR(i) = integral(@(x)integral_BR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauPFR(i) = integral(@(x)integral_PFR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauCSTR(i) = integral_CSTR(X(i),n,eps) / (k1*CA0^(n-1));
end

subplot(2,2,3);
loglog(1-X,TauBR, 1-X,TauPFR,'.', 1-X,TauCSTR);
title('eps = -0.5');
xlim([0.01 1]);
xlabel('1-X');
ylabel('residence time [s]');
legend('Batch','Plug Flow', 'CSTR');

% ------------------------------------------------------------------------%
% 4. Comparison PFR/CSTR
% ------------------------------------------------------------------------%
X = linspace(1e-2,0.999);

eps = 0;
for i=1:100
TauPFR_0(i) = integral(@(x)integral_PFR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauCSTR_0(i) = integral_CSTR(X(i),n,eps) / (k1*CA0^(n-1));
end

eps = 1;
for i=1:100
TauPFR_1(i) = integral(@(x)integral_PFR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauCSTR_1(i) = integral_CSTR(X(i),n,eps) / (k1*CA0^(n-1));
end

eps = -0.5;
for i=1:100
TauPFR_m05(i) = integral(@(x)integral_PFR(x,n,eps),0,X(i)) / (k1*CA0^(n-1));
TauCSTR_m05(i) = integral_CSTR(X(i),n,eps) / (k1*CA0^(n-1));
end

subplot(2,2,4);
loglog(1-X,TauCSTR_0./TauPFR_0,1-X,TauCSTR_1./TauPFR_1,1-X,TauCSTR_m05./TauPFR_m05);
title('CSTR vs PFR');
xlim([0.01 1]);
ylim([1 100]);
xlabel('1-X');
ylabel('residence time [s]');
legend('eps=0','eps=1','eps=-0.5');

% ------------------------------------------------------------------------%
% 5. Levenspiel's plots
% ------------------------------------------------------------------------%
figure;
hold on;
eps = 0;

n=1;
RA = -k1*CA0^n*(1-X).^n./(1+eps*X).^n;
loglog(X,-1./RA,'r', 'LineWidth',2);
n=2;
RA = -k1*CA0^n*(1-X).^n./(1+eps*X).^n;
loglog(X,-1./RA,'b','LineWidth',2);
n=0.5
RA = -k1*CA0^n*(1-X).^n./(1+eps*X).^n;
loglog(X,-1./RA,'g','LineWidth',2);

title('Levenspiel''s Plot');
xlim([0.01 1]);
ylim([1 100]);
xlabel('X (Conversion)');
ylabel('-1/R_A [m3.s/kmol]');
legend('n=1', 'n=2', 'n=0.5');