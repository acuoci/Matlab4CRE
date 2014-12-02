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

% The residence times for a single, reversible reaction
% in different reactors: batch (BR), plug flow (FPR), and CSTR:
%
% A <-> rB
% r = k1*(CA-CB/KC)
%
% The residence times are calculated through numerical integration of
% design equations.  

close all;
clear all;

% ------------------------------------------------------------------------%
% Kinetic constant [kmol, m3, s] (the units depend on the reaction order n)
% ------------------------------------------------------------------------%
k1 = 1;

% ------------------------------------------------------------------------%
% Initial concentration [kmol/m3/s]
% ------------------------------------------------------------------------%
CA0 = 1;

% ------------------------------------------------------------------------%
% Kinetic constant
% ------------------------------------------------------------------------%
prompt = 'What is the equilibrium constant? ';
KC = input(prompt);

% ------------------------------------------------------------------------%
% Relative initial amount of B (TetaB0)
% ------------------------------------------------------------------------%
prompt = 'What is the relative initial amount of B (TetaB0)? ';
tetaB0 = input(prompt);

% ------------------------------------------------------------------------%
% Conversion at equilibrium
% ------------------------------------------------------------------------%
eqX = @(KC,tetaB0,r) (KC-tetaB0)/(r+KC);

% ------------------------------------------------------------------------%
% Coefficient
% ------------------------------------------------------------------------%
Gamma = @(k1,KC,r) 1/k1*KC/(r+KC);

% ------------------------------------------------------------------------%
% Design equations for ideal reactors (functions to be integrated)
% ------------------------------------------------------------------------%
integrand_BR = @(x,  gamma,Xe,eps)  gamma* (Xe./(Xe-x));
integrand_PFR = @(x, gamma,Xe,eps)  gamma* (1+eps*x)./(Xe-x);
integral_CSTR = @(x, gamma,Xe,eps)  gamma* (1+eps*x)*x./(Xe-x);

% ------------------------------------------------------------------------%
% Conversion at equilibrium 
% ------------------------------------------------------------------------%


% ------------------------------------------------------------------------%
% Create plots
% ------------------------------------------------------------------------%
figure;                   % create new figure
xlimInf = 0.9;

% ------------------------------------------------------------------------%
% 1. Case 1: r=1 (no changes in number of moles)
% ------------------------------------------------------------------------%
r = 1;
gamma = Gamma(k1,KC,r);
Xe = eqX(KC,tetaB0,r);
eps = (r-1)/(1+tetaB0);
X = linspace(0,Xe-1e-2);  % create the x axis for plotting profiles
for i=1:100
TauBR(i) = integral(@(x)integrand_BR(x,gamma,Xe,eps),0,X(i));
TauPFR(i) = integral(@(x)integrand_PFR(x,gamma,Xe,eps),0,X(i));
TauCSTR(i) = integral_CSTR(X(i),gamma,Xe,eps);
end

subplot(2,2,1);
loglog(1-X,TauBR, 1-X,TauPFR,'.', 1-X,TauCSTR);
title('r = 1');
xlim([1-xlimInf, 1]);
xlabel('1-X');
ylabel('residence time [s]');
legend('Batch','Plug Flow', 'CSTR');

% ------------------------------------------------------------------------%
% 2. Case 2: r = 2 ( A<->2B)
% ------------------------------------------------------------------------%
r = 2;
gamma = Gamma(k1,KC,r);
Xe = eqX(KC,tetaB0,r);
eps = (r-1)/(1+tetaB0);
X = linspace(0,Xe-1e-2);  % create the x axis for plotting profiles
for i=1:100
TauBR(i) = integral(@(x)integrand_BR(x,gamma,Xe,eps),0,X(i));
TauPFR(i) = integral(@(x)integrand_PFR(x,gamma,Xe,eps),0,X(i));
TauCSTR(i) = integral_CSTR(X(i),gamma,Xe,eps);
end

subplot(2,2,2);
loglog(1-X,TauBR, 1-X,TauPFR,'.', 1-X,TauCSTR);
title('r = 2');
xlim([1-xlimInf, 1]);
xlabel('1-X');
ylabel('residence time [s]');
legend('Batch','Plug Flow', 'CSTR');

% ------------------------------------------------------------------------%
% 3. Case 3: r = 0.5
% ------------------------------------------------------------------------%
r = 0.5;
gamma = Gamma(k1,KC,r);
Xe = eqX(KC,tetaB0,r);
eps = (r-1)/(1+tetaB0);
X = linspace(0,Xe-1e-2);  % create the x axis for plotting profiles
for i=1:100
TauBR(i) = integral(@(x)integrand_BR(x,gamma,Xe,eps),0,X(i));
TauPFR(i) = integral(@(x)integrand_PFR(x,gamma,Xe,eps),0,X(i));
TauCSTR(i) = integral_CSTR(X(i),gamma,Xe,eps);
end

subplot(2,2,3);
loglog(1-X,TauBR, 1-X,TauPFR,'.', 1-X,TauCSTR);
title('r = 0.5');
xlim([1-xlimInf, 1]);
xlabel('1-X');
ylabel('residence time [s]');
legend('Batch','Plug Flow', 'CSTR');

% ------------------------------------------------------------------------%
% 4. Comparison PFR/CSTR
% ------------------------------------------------------------------------%
X = linspace(1e-2,0.99);

r = 1;
gamma = Gamma(k1,KC,r);
Xe = eqX(KC,tetaB0,r);
eps = (r-1)/(1+tetaB0);
for i=1:100
    TauPFR_0(i) = integral(@(x)integrand_PFR(x,gamma,Xe,eps),0,X(i));
    TauCSTR_0(i) = integral_CSTR(X(i),gamma,Xe,eps);
end

r = 2;
gamma = Gamma(k1,KC,r);
Xe = eqX(KC,tetaB0,r);
eps = (r-1)/(1+tetaB0);
for i=1:100
	TauPFR_1(i) = integral(@(x)integrand_PFR(x,gamma,Xe,eps),0,X(i));
    TauCSTR_1(i) = integral_CSTR(X(i),gamma,Xe,eps);
end

r = 0.5;
gamma = Gamma(k1,KC,r);
Xe = eqX(KC,tetaB0,r);
eps = (r-1)/(1+tetaB0);
for i=1:100
    TauPFR_m05(i) = integral(@(x)integrand_PFR(x,gamma,Xe,eps),0,X(i));
    TauCSTR_m05(i) = integral_CSTR(X(i),gamma,Xe,eps);
end

subplot(2,2,4);
loglog(1-X,TauCSTR_0./TauPFR_0,1-X,TauCSTR_1./TauPFR_1,1-X,TauCSTR_m05./TauPFR_m05);
title('CSTR vs PFR');
xlim([1-xlimInf, 1]);
ylim([1 10]);
xlabel('1-X');
ylabel('TauCSTR/TauPFR');
legend('r=1','r=2','r=0.5');

% ------------------------------------------------------------------------%
% 5. Levenspiel's plots
% ------------------------------------------------------------------------%
figure;
hold all;

r = 1;
eps = (r-1)/(1+tetaB0);

KC=1;
Xe = eqX(KC,tetaB0,r);
X = linspace(1e-2,Xe-1e-3);
RA = -k1*CA0./(1+eps*X) * (1+KC)/KC .* (Xe-X);
loglog(X,abs(-1./RA),'r');

KC=10;
Xe = eqX(KC,tetaB0,r);
X = linspace(1e-2,Xe-1e-3);
RA = -k1*CA0./(1+eps*X) * (1+KC)/KC .* (Xe-X);
loglog(X,-1./RA,'b');

KC=0.1;
Xe = eqX(KC,tetaB0,r);
X = linspace(1e-2,Xe-1e-3);
RA = -k1*CA0./(1+eps*X) * (1+KC)/KC .* (Xe-X);
loglog(X,-1./RA,'g');

title('Levenspiel''s Plot');
xlim([1e-2, 1]);
ylim([1 100]);
xlabel('X');
ylabel('-1/R_A');
legend('KC=1', 'KC=10', 'KC=0.1');