%% Main Script

% Gather User Input for Motor Components
[motorComponents, motorOperationalParams] = MotorThermalAnalysis.getUserInput();

% Extract operational parameters
V = motorOperationalParams.Voltage;
Kv = motorOperationalParams.Kv;
R_phase = motorOperationalParams.PhaseResistance;

% Define low current periods
lowCurrentStarts = [0, 300, 1000];
lowCurrentEnds = [120, 400, 1300];

% Simulation Parameters
total_time = 20000; % Total simulation time in seconds
num_steps = 60000;  % Number of simulation steps
time_step = total_time / num_steps; % Time step
T_environment = 24; % Ambient temperature (C)
time = linspace(0, total_time, num_steps); % Time array

% Initialize temperature arrays
T_rotor = zeros(1, num_steps);
T_stator = zeros(1, num_steps);
T_stator_simplified = zeros(1, num_steps); % For simplified stator model
T_initial_rotor = 25; % Initial temperature in °C
T_initial_stator = 25; % Initial temperature in °C
T_rotor(1) = T_initial_rotor;
T_stator(1) = T_initial_stator;
T_stator_simplified(1) = T_initial_stator; % For simplified stator model

% Simplified Thermal Parameters
C = motorComponents.stator.C_HeatCapacity;  % Total thermal capacitance (J/K)
R_thermal = 4.9163;  % Thermal resistance (K/W)
surfaceArea = pi * (motorComponents.stator.OuterDiameter / 1000) * (motorComponents.stator.Length / 1000); % Surface area of the stator in m^2

%% Main Loop for Stator and Rotor Temperature Analysis

for i = 2:num_steps
    % Dynamic current value based on gait profile
    current = MotorThermalAnalysis.currentProfileForGait(time(i), total_time, lowCurrentStarts, lowCurrentEnds);

    % Example relationship for motor speed
    someFactor = 250; % Adjust this factor based on your motor's characteristics
    motorSpeed = someFactor * current; % Calculate dynamic motor speed

    % Calculate Power Loss in Stator
    powerLoss = current^2 * R_phase;

    % Calculate Heat Transfer Calculations for Stator
    Q_generated_stator = powerLoss; % Heat generated in stator
    dTdt_stator = (Q_generated_stator - (T_stator(i-1) - T_environment) / R_thermal) / C;
    T_stator(i) = T_stator(i-1) + dTdt_stator * time_step;

    % Calculate convective heat transfer coefficient for simplified stator model
    h = MotorThermalAnalysis.calcConvectiveHeatTransferCoeff(motorSpeed, motorComponents.stator.OuterDiameter, T_stator_simplified(i-1), T_environment);
    R_conv = MotorThermalAnalysis.calcConvectiveResistance(h, surfaceArea);

    % Update Simplified Stator Temperature
    Q_lost_conv = (T_stator_simplified(i-1) - T_environment) / R_conv;
    dTdt_stator_simplified = (Q_generated_stator - Q_lost_conv) / C;
    T_stator_simplified(i) = T_stator_simplified(i-1) + dTdt_stator_simplified * time_step;

    % Calculate Heat Transfer from Stator to Rotor
    Q_transfer = MotorThermalAnalysis.estimateHeatTransfer(motorComponents.axle, motorComponents.stator, T_stator(i), T_rotor(i-1));

    % Update Rotor Temperature based on Heat Transfer
    dTdt_rotor = Q_transfer / motorComponents.rotor.C_ThermalCapacitance;
    T_rotor(i) = T_rotor(i-1) + dTdt_rotor * time_step;
end

%% Plotting Results
time = 0:time_step:(total_time-time_step);
plot(time, T_stator, 'b', 'LineWidth', 1);
hold on;  % Hold on to the current plot
plot(time, T_rotor, 'r', 'LineWidth', 1);
plot(time, T_stator_simplified, 'g', 'LineWidth', 1);
hold off; % Release the plot hold

xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Stator and Rotor Thermal Simulation');
legend('Stator - Complex Model', 'Rotor - Complex Model', 'Stator - Simplified Model');
grid on;
