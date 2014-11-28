clear all;

% ------------------------------------------------------------------------%
% Global variables 
% ------------------------------------------------------------------------%
global kinetics;            % kinetic mechanism
global T0;                  % inlet temperature [K]
global P0;                  % inlet pressure [Pa]
global omega0;              % inlet mass fractions [-]
global v0;                  % inlet velocity [m/s]
global D;                   % internal diameter
global gx;                  % gravity [m/s2]
global U;
global Te;

% ------------------------------------------------------------------------%
% 1. Load the kinetic mechanism 
% ------------------------------------------------------------------------%
kinetics = GasMechanism_AceticAnhydride;

% ------------------------------------------------------------------------%
% 2. Define input data
% ------------------------------------------------------------------------%
T0 = 1000.;                     % inlet temperature [K]
P0 = 101325;                    % inlet pressure [Pa]

omega0 = zeros(1,kinetics.ns);  % inlet mass fractions [-]
omega0(1) = 1.;

D = 2.66e-2;                    % internal diameter [m]
v0 = 10;                        % inlet velocity [m/s]
x0 = 0;                         % initial axial coordinate [m]
xF = 10;                        % final axial coordinate [m]

U  = 30;                        % global heat exchange coefficient [W/m2/K]
Te = 1100;                      % external temperature [K]

gx = 9.81;                      % gravity (alog the axis) [m/s2]
                                % horizontal: gx = 0
                                % vertical, going up: gx = -9.81 m/s2
                                % vertical, going down: gx =  9.81 m/s2
                                
                                
% ------------------------------------------------------------------------%
% 3. Simulation
% ------------------------------------------------------------------------%
clear screen;
disp( sprintf('1. Isothermal, without pressure drop'));
disp( sprintf('2. Isothermal, with pressure drop'));
disp( sprintf('3. Non-Isothermal, without pressure drop'));
disp( sprintf(' '));
prompt = '   What type of simulation? ';
simulation = input(prompt);


if (simulation == 1)
    
    % Initial conditions
    y0 = [omega0, T0, P0];
    
    % ODE Solver
    [x,y] = ode15s(@odeTubularIsothermal,[x0 xF], y0);
    
elseif (simulation == 2)

    % Initial conditions
    y0 = [omega0, T0, P0];
    
    % ODE Solver
    [x,y] = ode15s(@odeTubularIsothermalWithPressureDrop,[x0 xF], y0);
    
elseif (simulation == 3)
    
    % Initial conditions
    y0 = [omega0, T0, P0];
    
    % ODE Solver
    [x,y] = ode15s(@odeTubularHeatExchange,[x0 xF], y0);
    
end

% ------------------------------------------------------------------------%
% 3. Post processing
% ------------------------------------------------------------------------%

figure; % create new figure
hold all;

% Mass fraction profiles
for i=1:kinetics.ns
    subplot(2,2,1);
    hold all;
    plot (x, y(:,i),'LineWidth',2);
end
title('Mass fractions of species along the reactor');
xlabel('axial length [m]'); 
ylabel('mass fractions [-]'); 
legend(kinetics.species);
axis([-inf,inf,0,inf]);

% Conversion profile
subplot(2,2,2);
plot (x, 1-y(:,1),'LineWidth',2);
title('Conversion profile along the reactor');
xlabel('axial length [m]'); 
ylabel('Conversion');

% Temperature profile
subplot(2,2,3);
plot (x, y(:,kinetics.ns+1),'LineWidth',2);
title('Temperature profile along the reactor');
xlabel('axial length [m]'); 
ylabel('Temperature [K]');

% Pressure profile
subplot(2,2,4);
plot (x, y(:,kinetics.ns+2)/1e5,'LineWidth',2);
title('Pressure profile along the reactor');
xlabel('axial length [m]'); 
ylabel('Pressure [bar]');
