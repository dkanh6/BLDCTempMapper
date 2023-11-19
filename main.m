%% Main Script

% Gather User Input for Motor Components
[motorComponents, motorOperationalParams] = MotorThermalAnalysis.getUserInput();

% Extract operational parameters
V = motorOperationalParams.Voltage;
I = motorOperationalParams.Current;
Kv = motorOperationalParams.Kv;
R_phase = motorOperationalParams.PhaseResistance;

% Estimate Motor Parameters (Speed and Power Loss)
[motorSpeed, powerLoss] = MotorThermalAnalysis.estimateMotorParameters(V, I, Kv, R_phase);

% Display the estimated motor speed and power loss
disp(['<strong>Motor Speed:</strong> ', num2str(motorSpeed), ' RPM']);
disp(['<strong>Power Loss:</strong> ', num2str(powerLoss), ' Watts']);

%% Calculate Thermal Resistances
% Calculate Conductive Equivalent Resistances
R_thermal_rotor = MotorThermalAnalysis.calcThermalResistance(motorComponents.rotor);
R_thermal_stator = MotorThermalAnalysis.calcThermalResistance(motorComponents.stator);

% Calculate Convective Equivalent Resistances
statorDiameter = motorComponents.stator.OuterDiameter;
rotorDiameter = motorComponents.rotor.OuterDiameter;
h_stator = MotorThermalAnalysis.calcConvectiveHeatTransferCoeff(motorSpeed, statorDiameter);
h_rotor = MotorThermalAnalysis.calcConvectiveHeatTransferCoeff(motorSpeed, rotorDiameter);
surface_area_stator = 2 * pi * (motorComponents.stator.OuterDiameter / 1000) * (motorComponents.stator.Length/1000);
surface_area_rotor = 2 * pi * (motorComponents.rotor.OuterDiameter / 1000) * (motorComponents.rotor.Length/1000);
R_conv_stator_rotor = MotorThermalAnalysis.calcConvectiveResistance(h_stator, surface_area_stator);
R_conv1_rotor_env = MotorThermalAnalysis.calcConvectiveResistance(h_rotor, surface_area_rotor);

% Assume initial temperatures
T_initial_rotor = 25; % Initial temperature in °C
T_initial_stator = 25; % Initial temperature in °C

% Simulation Parameters
total_time = 3000; % Total simulation time in seconds
num_steps = 300;  % Number of simulation steps

% Time array for simulation
time = linspace(0, total_time, num_steps);

%% Thermal Model Simulation
for i = 2:num_steps
    % Heat generated in the stator (Assume all power loss converts to heat)
    Q_generated = powerLoss; % Heat generated in stator (Watts)

    % Heat transfer from stator to rotor through axle (Assume simple linear conduction)
    delta_T_axle = T_stator(i-1) - T_rotor(i-1); % Temperature difference
    Q_conduction = delta_T_axle / R_thermal_rotor; % Conduction through rotor (Watts)

    % Heat dissipation from stator and rotor to environment
    Q_dissipation_stator = (T_stator(i-1) - T_environment) / R_conv_stator_rotor; % Stator (Watts)
    Q_dissipation_rotor = (T_rotor(i-1) - T_environment) / R_conv1_rotor_env; % Rotor (Watts)

    % Update temperatures
    T_stator(i) = T_stator(i-1) + time_step * (Q_generated - Q_conduction - Q_dissipation_stator) / (motorComponents.stator.C_HeatCapacity);
    T_rotor(i) = T_rotor(i-1) + time_step * (Q_conduction - Q_dissipation_rotor) / (motorComponents.rotor.C_HeatCapacity);

    % Check for Tmax
    if T_stator(i) > motorOperationalParams.Tmax || T_rotor(i) > motorOperationalParams.Tmax
        warning('Maximum temperature exceeded at time %f seconds', time(i));
        break;
    end
end

%% Plotting Results
plot(time, T_rotor, 'r', time, T_stator, 'b');
legend('Rotor Temperature', 'Stator Temperature');
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Motor Thermal Simulation');
grid on;
