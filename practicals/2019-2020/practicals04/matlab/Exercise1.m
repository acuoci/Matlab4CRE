close all; clear variables;

%% Input data

% Operating conditions
P0   = 6*101325; % pressure [Pa]
T0   = 373;      % temperature [K]
km   = 0.2;      % mass transfer coefficient [1/min]
FA0  = 15;       % inlet flow rate of species A [mol/min]
FB0  =  0;       % inlet flow rate of species B [mol/min]
FC0  =  0;       % inlet flow rate of species C [mol/min]
Vtot = 1000;     % final volume [dm3]

% Kinetic and thermodynamic data
k = 0.7;   % forward kinetic constant [1/min]                  
KC = 0.05; % equilibrium constant[mol/dm3]  

% Initial conditions
F0 = [FA0, FB0, FC0];


%% Case 1: pure plug flow reactor

% ODE solution
[V,F] = ode15s(@plugflow,[0 Vtot], F0, [], k, KC, T0, P0);

% Plotting solution
figure;
subplot(2,3,1);
hold all;
plot(V,F(:,1), 'LineWidth', 2);
plot(V,F(:,2), 'LineWidth', 2);
plot(V,F(:,3), 'LineWidth', 2, 'LineStyle', '--');
title ('Pure plug flow reactor'); legend('A','B','C');
xlabel('volume [dm^3]'); ylabel('flow rates [mol/min]');

subplot(2,3,2);
hold all;
Ctot = P0/8314/T0;
Ftot = F(:,1)+F(:,2)+F(:,3);
plot(V,Ctot*F(:,1)./Ftot, 'LineWidth', 2);
plot(V,Ctot*F(:,2)./Ftot, 'LineWidth', 2);
plot(V,Ctot*F(:,3)./Ftot, 'LineWidth', 2, 'LineStyle', '--');
title ('Pure plug flow reactor');   legend('A','B','C');
xlabel('volume [dm^3]'); ylabel('concentrations [mol/dm3]');

subplot(2,3,3);
hold all;
[hAx,hLine1,hLine2] = plotyy(V,Ftot, V, (FA0-F(:,1))/FA0);
title ('Pure plug flow reactor'); legend('Total molar flow rate', 'Conversion');
xlabel('volume [dm^3]');
set(hLine1, 'LineWidth', 2); set(hLine2, 'LineWidth', 2);
ylabel(hAx(1),'total molar flow rate [mol/min]') % left y-axis
ylabel(hAx(2),'conversion [-]')                  % right y-axis


%% Case 2: membrane reactor (only B component)

% ODE solution
[V,F] = ode15s(@membranereactor_only_B,[0 Vtot], F0, [], k,KC,T0,P0,km);

% Plotting solution
subplot(2,3,4);
hold all;
plot(V,F(:,1), 'LineWidth', 2);
plot(V,F(:,2), 'LineWidth', 2);
plot(V,F(:,3), 'LineWidth', 2, 'LineStyle', '--');
title ('Membrane reactor (only B)'); legend('A','B','C');
xlabel('volume [dm^3]'); ylabel('flow rates [mol/min]');

subplot(2,3,5);
hold all;
Ctot = P0/8314/T0;
Ftot = F(:,1)+F(:,2)+F(:,3);
plot(V,Ctot*F(:,1)./Ftot, 'LineWidth', 2);
plot(V,Ctot*F(:,2)./Ftot, 'LineWidth', 2);
plot(V,Ctot*F(:,3)./Ftot, 'LineWidth', 2, 'LineStyle', '--');
title ('Membrane reactor (only B)'); legend('A','B','C');
xlabel('volume [dm^3]'); ylabel('concentrations [mol/dm^3]');

subplot(2,3,6);
hold all;
[hAx,hLine1,hLine2] = plotyy(V,Ftot, V, (FA0-F(:,1))/FA0);
legend('Total molar flow rate', 'Conversion'); title ('Membrane reactor (only B)');
xlabel('volume [dm^3]');
set(hLine1, 'LineWidth', 2); set(hLine2, 'LineWidth', 2);
ylabel(hAx(1),'total molar flow rate [mol/min]') % left y-axis
ylabel(hAx(2),'conversion [-]')                  % right y-axis



%% ODE system 1: pure plug-flow reactor
function dF = plugflow(t,F, k,KC,T0,P0)
        
    FA = F(1);                % molar flow rate of species A [mol/min]
    FB = F(2);                % molar flow rate of species A [mol/min]
    FC = F(3);                % molar flow rate of species A [mol/min]

    Ftot = FA+FB+FC;          % total molar flow rate [mol/min]    
    Ctot = P0/8314/T0;        % total concentration [mol/dm3]
    CA = Ctot*FA/Ftot;        % concentration of A [mol/dm3]
    CB = Ctot*FB/Ftot;        % concentration of A [mol/dm3]
    CC = Ctot*FC/Ftot;        % concentration of A [mol/dm3]
    
    r  = k*(CA-CB*CC/KC);     % reaction rate [mol/dm3/min]

    dFA = -r;                 % dFA_over_dV [mol/dm3/min]
    dFB =  r;                 % dFB_over_dV [mol/dm3/min]
    dFC =  r;                 % dFC_over_dV [mol/dm3/min]
    
    dF = [dFA dFB dFC]';
    
end


%% ODE system 2: membrane reactor (only B component)
function dF = membranereactor_only_B(t,F, k,KC,T0,P0,km)
    
    FA = F(1);                % molar flow rate of species A [mol/min]
    FB = F(2);                % molar flow rate of species B [mol/min]
    FC = F(3);                % molar flow rate of species C [mol/min]
    
    Ftot = FA+FB+FC;          % total molar flow rate [mol/min]    
    Ctot = P0/8314/T0;        % total concentration [mol/dm3]
    CA = Ctot*FA/Ftot;        % concentration of A [mol/dm3]
    CB = Ctot*FB/Ftot;        % concentration of A [mol/dm3]
    CC = Ctot*FC/Ftot;        % concentration of A [mol/dm3]
    
    r  = k*(CA-CB*CC/KC);     % reaction rate [mol/dm3/min]
    
    dFA = -r;                 % dFA_over_dV [mol/dm3/min]
    dFB =  r - km*CB;         % dFB_over_dV [mol/dm3/min]
    dFC =  r;                 % dFC_over_dV [mol/dm3/min]
    
    dF = [dFA dFB dFC]';
end
