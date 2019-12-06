clear variables; close all;

global A1 A2;
global E1 E2;
global DH1 DH2;
global tau Cp;
global V0 Q0 CA0 T0; 
global U A Text;

%% Input data

Cp = 200;       % specific heat [J/mol/K]
CA0 = 0.3;      % inlet concentration [mol/l]
T0 = 283;       % inlet temperature [K]
V0 = 10;        % volume [l]
Q0 = 1000;      % volumetric flow rate [l/min]

k1T1 = 3.30;    % kinetic constant k1@T1 [1/min]
k2T2 = 4.58;    % kinetic constant k2@T2 [1/min]
T1 = 300;       % temperature for k1T1 [K]
T2 = 500;       % temperature for k2T2 [K]
E1 = 9900;      % activation energy [cal/mol]
E2 = 27000;     % activation energy [cal/mol]
DH1 = -55000;   % reaction enthalpy [J/mol]
DH2 = -71500;   % reaction enthalpy [J/mol]

CpIn = Cp;      % specific heat [J/mol/K]
U = 4000;       % heat exchange coefficient [J/min/K/m2]
A = 10;         % exchange area [m2]
Text = 57+273;  % external fluid temperature [K]


%% Pre-processing of input data

% Calculation of reaction frequency factors
A1 = k1T1/exp(-E1/1.987/T1);    % [1/min]
A2 = k2T2/exp(-E2/1.987/T2);    % [1/min]

% Calculation of additional quantities
tau = V0/Q0;                        % residence time [min]
FA0 = CA0*Q0;                       % molar flow rate [mol/min]
kappa = U*A/FA0/CpIn;               % heat exchange
Tstar = (T0+kappa*Text)/(1+kappa);  % heat exchange

% Definition of temperature range
T = 300:1:700;  % [K]

% Kinetic constants as a function of T
k1 = A1*exp(-E1/1.987./T);  % [1/min]
k2 = A2*exp(-E2/1.987./T);  % [1/min]

% Concentrations as a function of T
CA = CA0./(1+k1*tau);           % [mol/l]
CB = CA.*k1*tau./(1+k2*tau);    % [mol/l]

% Reaction rates as a function of T
r1 = k1.*CA;    % [mol/l/min]
r2 = k2.*CB;    % [mol/l/min]


%% Calculation and plotting of R and G curves

% Removed and generated heats
R = CpIn*(1+kappa)*(T-Tstar);          % [J/mol]
G = -DH1*r1*V0/FA0 -DH2*r2*V0/FA0;     % [J/mol]

% Plot the lines to find the solution
figure;
plot(T,R,T,G);
legend('R', 'G', 'Location', 'NorthWest'); 
axis([300 700 0 inf]);
xlabel('temperature [K]'); ylabel('R, G [J/mol]');
grid on;

%% Numerically find the intersections
delta = abs(G-R);
TSS = zeros(1,5);
for i = 1:5
    [val,index] = min(delta);
    delta(index) = max(delta);
    TSS(i) = T(index);
end
[minval,index] = min(delta);

% Compute concentrations at the different SS temperatures
k1SS = A1.*exp(-E1./1.987./TSS); % [l/min]
k2SS = A2.*exp(-E2./1.987./TSS);
CASS = CA0./(1+k1SS.*tau);
CBSS = CASS.*k1SS.*tau./(1+k2SS.*tau);
CCSS = CBSS.*k2SS.*tau;

figure; hold on; grid on;
plot(TSS,CASS,'or');
plot(TSS,CBSS,'ob');
plot(TSS,CCSS,'og');
hold off; xlabel('T [K]'); ylabel('C [mol/l]');
legend('CA','CB','CC','Location','SouthEast');
title('Concentrations at SS conditions');


%% Part 2: solve the system with fsolve
%T_guess = T0;
%TSS_fsolve = fsolve(@heat_balance, T_guess);


%% Part 3: solve the system as an unsteady CSTR
% the function is the same, but you have to write it 
% as a function of time
tf = 1; % assume the final time as 1 minute
[t,T_unsteady] = ode15s(@heat_balance_unsteady,[0 tf],[T0]);

figure; grid on;
plot(t,T_unsteady,'LineWidth',2)
title('Temperature Profile in the CSTR')
xlabel('time [min]'); ylabel('T [K]');

% Consider the fluctuation of T0 in time
[t_T0var,T_T0var] = ode15s(@heat_balance_unsteady_T0var,[0 tf],T0);
T0var = 283 + 20.*((0.5<t_T0var).*(t_T0var<0.55));

figure; grid on;
plot(t_T0var,T_T0var,t_T0var,T0var,'LineWidth',2)
title('Temperature Profile in the CSTR - inlet T fluctuation')
xlabel('time [min]'); ylabel('T [K]');
legend('T CSTR','T INLET','Location','NorthWest')


%% Non Linear Equation (Energy balance equation)
function Q = heat_balance(T)

    global A1 A2;
    global E1 E2;
    global DH1 DH2;
    global tau Cp;
    global V0 Q0 CA0 T0; 
    global U A Text;

    % Updating kinetic constants
    k1 = A1.*exp(-E1./1.987./T);    % [l/min]
    k2 = A2.*exp(-E2./1.987./T);    % [l/min]

    % Concentrations as a function of temperature
    tau = V0/Q0;                    % [min]
    CA = CA0./(1+k1.*tau);          % [mol/l]
    CB = CA.*k1.*tau./(1+k2.*tau);  % [mol/l]
    CC = CB.*k2.*tau;               % [mol/l]

    % Heat generated G and heat removed R
    r1 = k1.*CA;        % [mol/l/min]
    r2 = k2.*CB;        % [mol/l/min]
    FA0 = CA0*Q0;       % [mol/min]
    G = -V0.*(r1.*DH1+r2.*DH2)./FA0;     % [J/mol]
    R = -Cp.*(T0-T)-U.*A.*(Text-T)./FA0; % [J/mol]

    % Output function to be set to 0
    Q = G-R;
    
end


%% Unsteady balance equation
function dQ = heat_balance_unsteady(t,T)

    global A1 A2;
    global E1 E2;
    global DH1 DH2;
    global tau Cp;
    global V0 Q0 CA0 T0; 
    global U A Text;
    
    % Updating kinetic constants
    k1 = A1.*exp(-E1./1.987./T);    % [l/min]
    k2 = A2.*exp(-E2./1.987./T);    % [l/min]

    % Concentrations as a function of temperature
    tau = V0/Q0;                    % [min]
    CA = CA0./(1+k1.*tau);          % [mol/l]
    CB = CA.*k1.*tau./(1+k2.*tau);  % [mol/l]
    CC = CB.*k2.*tau;               % [mol/l]

    % Heat generated G and the heat removed R
    r1 = k1.*CA;                % [mol/l/min]
    r2 = k2.*CB;                % [mol/l/min]
    FA0 = CA0*Q0;               % [mol/min]
    G = -V0.*(r1.*DH1+r2.*DH2)./FA0;        % [J/mol]
    R = -Cp.*(T0-T)-U.*A.*(Text-T)./FA0;    % [J/mol]

    % Output function to be set to 0
    dQ = G-R;
    
end


%% Unsteady balance equation with disturbance
function dQ = heat_balance_unsteady_T0var(t,T)

    global A1 A2;
    global E1 E2;
    global DH1 DH2;
    global tau Cp;
    global V0 Q0 CA0 T0; 
    global U A Text;

    % Disturbance
    T0 = 283 + 20.*((0.5<t).*(t<0.55));

    % Updating the kinetic constants
    k1 = A1.*exp(-E1./1.987./T);    % [l/min]
    k2 = A2.*exp(-E2./1.987./T);    % [l/min]

    % Concentrations as a function of temperature
    tau = V0/Q0;                    % [min]
    CA = CA0./(1+k1.*tau);          % [mol/l]
    CB = CA.*k1.*tau./(1+k2.*tau);  % [mol/l]
    CC = CB.*k2.*tau;               % [mol/l]

    % Heat generated G and the heat removed R
    r1 = k1.*CA;                            % [mol/l/min]
    r2 = k2.*CB;                            % [mol/l/min]
    FA0 = CA0*Q0;                           % [mol/min]
    G = -V0.*(r1.*DH1+r2.*DH2)./FA0;        % [J/mol]
    R = -Cp.*(T0-T)-U.*A.*(Text-T)./FA0;    % [J/mol]

    % Output function to be set to 0
    dQ = G-R;
    
end
