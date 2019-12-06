clear variables; close all;

global FSin;
global kappa;
global Qin;

%% Input data

V0 = 100;           % initial volume [l]
Qin = 5;            % inlet volumetric flow rate [l/min]
kappa = 0.05;       % kinetic constant [l/mol/min]
CSin = 2.0;         % inlet concentration of S [mol/l]
CA0 = 10;           % initial concentration of A [mol/l]
NS0 = 0;            % initial number of moles of S [mol]
tf = 30;            % batch time [min]


%% Solution

NA0 = CA0*V0;       % initial number of moles of A [mol]
FSin = CSin*Qin;    % inlet molar flow rate for S [mol/min]

% ODE solution
Yin = [NA0 NS0 V0]';
[t, Y] = ode45(@SemiBatch, [0 tf], Yin);


%% Post-processing

NA = Y(:,1);    % moles of A [mol]
NS = Y(:,2);    % moles of S [mol]
V  = Y(:,3);    % volume [l]

CA = NA./V;     % concentration of A [mol/l]
CS = NS./V;     % concentration of S [mol/l]

% Plot moles
figure;
plot(t,NA, t,NS); legend('NA', 'NS');
xlabel('time [s]'); ylabel('number of moles [mol]');

% Plot concentrations
figure;
plot(t,CA, t,CS); legend('CA', 'CS');
xlabel('time [s]'); ylabel('concentrations [mol/l]');


%% Semibatch ODE function
function dY = SemiBatch(t,Y)

    global FSin;
    global kappa;
    global Qin;

    NA = Y(1);
    NS = Y(2);
    V  = Y(3);
    
    CA = NA/V;
    CS = NS/V;
    
    r = kappa*CA*CS;
    dNA = -r*V;
    dNS = FSin;
    dV  = Qin;
    
    dY = [dNA dNS dV]';
    
end
