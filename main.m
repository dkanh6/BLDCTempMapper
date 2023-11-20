%% Main Script

% Gather User Input for Motor Components
[motorComponents, motorOperationalParams] = MotorThermalAnalysis.getUserInput();

% Extract operational parameters
V = motorOperationalParams.Voltage;
Kv = motorOperationalParams.Kv;
R_phase = motorOperationalParams.PhaseResistance;
baseSpeed = Kv * V;

% Define low current periods
lowCurrentStarts = [0, 300, 1020];
lowCurrentEnds = [120, 420, 1320];

% Simulation Parameters
total_time = 1000; % Total simulation time in seconds
num_steps = 1000;  % Number of simulation steps
time_step = total_time / num_steps; % Time step
T_environment = 24; % Ambient temperature (C)
time = linspace(0, total_time, num_steps); % Time array

% Initialize temperature arrays
T_rotor = zeros(1, num_steps);
T_stator = zeros(1, num_steps);
T_initial_rotor = 25; % Initial temperature in °C
T_initial_stator = 25; % Initial temperature in °C
T_rotor(1) = T_initial_rotor;
T_stator(1) = T_initial_stator;

% Simplified Thermal Parameters (Based on the simplified model)
C = motorComponents.stator.C_HeatCapacity;  % Total thermal capacitance (J/K)
R_thermal = 4.9163;  % Thermal resistance (K/W)

%% Main Loop Simplified for Stator and Rotor Temperature Analysis

for i = 2:num_steps
    % Dynamic current value based on the profile
    current = MotorThermalAnalysis.currentProfile(time(i), lowCurrentStarts, lowCurrentEnds);

    % Calculate Power Loss in Stator (Simplified)
    powerLoss = current^2 * R_phase;

    % Calculate Heat Transfer Calculations for Stator (Simplified)
    Q_generated_stator = powerLoss; % Heat generated in stator
    dTdt_stator = (Q_generated_stator - (T_stator(i-1) - T_environment) / R_thermal) / C;
    T_stator(i) = T_stator(i-1) + dTdt_stator * time_step;

    % Calculate Heat Transfer from Stator to Rotor
    Q_transfer = MotorThermalAnalysis.estimateHeatTransfer(motorComponents.axle, motorComponents.stator, T_stator(i), T_rotor(i-1));

    % Update Rotor Temperature based on Heat Transfer
    dTdt_rotor = Q_transfer / motorComponents.rotor.C_ThermalCapacitance;
    T_rotor(i) = T_rotor(i-1) + dTdt_rotor * time_step;

    % Display Information for Debugging
    disp(['Time: ', num2str(time(i)), ' s, Stator Temperature: ', num2str(T_stator(i)), ' °C, Rotor Temperature: ', num2str(T_rotor(i)), ' °C']);

    % Check for Unrealistic Temperatures
    if T_stator(i) > 1000 || T_rotor(i) > 1000 % Use realistic thresholds for your motor
        disp(['Warning: Exceedingly high temperature at time ', num2str(time(i)), ' s. Check calculations and parameters.']);
        break; % Break out of the loop if temperatures are unrealistically high
    end
end 

%% Plotting Results
plot(time, T_stator, 'b', time, T_rotor, 'r');
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Stator and Rotor Thermal Simulation');
legend('Stator', 'Rotor');
grid on;

%% Auxiliary Functions



