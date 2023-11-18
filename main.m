%% Main Script

% In your main script or function
motorComponents = MotorThermalAnalysis.getUserInput();

% Now motorComponents has stator and rotor fields that you can use
statorResistance = MotorThermalAnalysis.calcThermalResistance(motorComponents.stator);
rotorResistance = MotorThermalAnalysis.calcThermalResistance(motorComponents.rotor);




% Set up Paramters
R_thermal_rotor = 0;  
C_thermal_rotor = 0;       % Total thermal capacitance of the rotor (J/K)
R_thermal_stator = 0;
C_thermal_stator = 0;                % Total thermal apacitance of the stator (J/K) r_thermal_rotor = 1;         % Thermal resistance of the rotor (K/W)% Thermal resistance of the stator (K/W)
T_ambient = 24;              % Ambient temperature (°C)
total_time = 3000;           % Total time for simulation (s)
num_steps = 9000;            % Number of steps in the simulation
R_phase = 0.055178;         % Resistance of one phase
K_transfer = 0;            % Heat transfer coefficient (W/K)



