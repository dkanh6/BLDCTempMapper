
classdef MotorThermalAnalysis
    methods (Static)
        function [cp_composite, k_composite] = calcCompositeProperties(copperFraction)
            % Define the specific heat capacities (J/g°C)
            cp_copper = 385; % Copper
            cp_steel = 466;  % Steel

            % Define the thermal conductivities (W/mK)
            k_copper = 401;    % Copper
            k_steel = 50;      % Steel

            % Calculate the fraction of steel based on the given copper fraction
            steelFraction = 1 - copperFraction;

            % Composite specific heat capacity (cp)
            cp_composite = (cp_copper * copperFraction) + (cp_steel * steelFraction);

            % Composite thermal conductivity (k)
            % Using a simple linear rule of mixture (not strictly series or
            % parallell construction in the stator)
            k_composite = (k_copper * copperFraction) + (k_steel * steelFraction);
        end

        function C = calcThermalCapacitance(geometry, cp_composite, componentType)
            % Extract geometric properties from structure
            outerDiameter = geometry.OuterDiameter / 1000; % Convert to meters
            innerDiameter = geometry.InnerDiameter / 1000; % Convert to meters
            length = geometry.Length / 1000; % Convert to meters

            % Calculate the volume of the hollow cylinder (assuming it's a hollow cylinder)
            volume = pi * (outerDiameter^2 - innerDiameter^2) / 4 * length;

            % Set density based on component type
            switch componentType
                case 'stator'
                    % Density for stator (half copper and half steel)
                    density = (8000 + 8960) / 2; % Example: average of steel and copper density
                case 'rotor'
                    % Density for rotor (aluminum)
                    density = 2700; % Example: density of aluminum
                case 'axle'
                    % Density for axle (steel)
                    density = 8000; % Example: density of steel
                otherwise
                    error('Invalid component type');
            end

            % Calculate mass
            mass = density * volume;

            % Calculate thermal capacitance
            C = cp_composite * mass;
        end

        function [motorComponents, motorOperationalParams] = getUserInput()

            defaultTmax = 60;
            % Display default values
            disp('<strong>Default Motor Component and Operational Parameters:</strong>');
            disp('---------------------------------------------------');
            disp('<strong>Stator:</strong> Outer Diameter = 27.5 mm, Inner Diameter = 5 mm, Length = 26 mm, Thermal Conductivity = 225.5 W/mK, Heat Capacity = 0.4255 J/kgK');
            disp('<strong>Rotor:</strong> Outer Diameter = 35 mm, Inner Diameter = 30 mm, Length = 36 mm, Thermal Conductivity = 167 W/mK, Heat Capacity = 0.89 J/kgK');
            disp('<strong>Axle:</strong> Outer Diameter = 5 mm, Length = 55 mm, Thermal Conductivity = 50 W/mK, Heat Capacity = 0.466 J/kgK');
            disp('<strong>Operational Parameters:</strong> Voltage = 14.7 V, Kv = 910 RPM/V, ax Motor Temp = 60C Phase Resistance = 0.055178 Ohms');
            disp('---------------------------------------------------');
            % Ask user to choose between default values or custom values
            choice = MotorThermalAnalysis.getYorNInput('Do you want to use default values? (Y/N): ');

            if upper(choice) == 'Y'

                % Use default values
                % Example default values
                % cp copper = 0.385 J/gC
                % cp Aluminum = 0.89 J/gC
                % cp Steel = 0.466 J/gC
                [cp_stator, k_stator] = MotorThermalAnalysis.calcCompositeProperties(0.5); % 50% copper
                stator = struct('OuterDiameter', 27.5, 'InnerDiameter', 5, 'Length', 26, 'K_Conductivity', k_stator, 'C_HeatCapacity', cp_stator); % Assuming 50:50 ratio of copper to steel to calculate the thermal capacitance and conductivity
                rotor = struct('OuterDiameter', 35, 'InnerDiameter', 30, 'Length', 36, 'K_Conductivity', 167, 'C_HeatCapacity', 890);
                axle = struct('OuterDiameter', 5, 'InnerDiameter', 0.1, 'Length', 55, 'K_Conductivity', 50, 'C_HeatCapacity', 466);
                motorComponents = struct('stator', stator, 'rotor', rotor, 'axle', axle);
                motorOperationalParams = struct('Voltage', 14.7, 'Current', 10, 'Kv', 910, 'PhaseResistance', 0.055178, 'Tmax',defaultTmax);
                motorComponents.stator.C_ThermalCapacitance = MotorThermalAnalysis.calcThermalCapacitance(motorComponents.stator, motorComponents.stator.C_HeatCapacity, 'stator');
                motorComponents.rotor.C_ThermalCapacitance = MotorThermalAnalysis.calcThermalCapacitance(motorComponents.rotor, motorComponents.rotor.C_HeatCapacity, 'rotor');
                motorComponents.axle.C_ThermalCapacitance = MotorThermalAnalysis.calcThermalCapacitance(motorComponents.axle, motorComponents.axle.C_HeatCapacity, 'axle');
            else
                % Prompt for custom values
                disp('<strong>==============================</strong>');
                disp('<strong>CUSTOM MOTOR COMPONENT INPUT</strong>');
                disp('<strong>==============================</strong>');

                % Stator Information
                disp('<strong>STATOR INFORMATION:</strong>');
                disp('-------------------');
                stator.OuterDiameter = input('Enter Stator Outer Diameter (mm): ');
                stator.InnerDiameter = input('Enter Stator Inner Diameter (mm): ');
                stator.Length = input('Enter Stator Length (mm): ');
                stator.K_Conductivity = input('Enter Stator Thermal Conductivity (W/mK): ');
                stator.C_HeatCapacity = input('Enter Stator Specific Heat Capacity (J/kgK): ');
                disp(' ');

                % Rotor Information
                disp('<strong>ROTOR INFORMATION:</strong>');
                disp('------------------');
                rotor.OuterDiameter = input('Enter Rotor Outer Diameter (mm): ');
                rotor.InnerDiameter = input('Enter Rotor Inner Diameter (mm): ');
                rotor.Length = input('Enter Rotor Length (mm): ');
                rotor.K_Conductivity = input('Enter Rotor Thermal Conductivity (W/mK): ');
                rotor.C_HeatCapacity = input('Enter Rotor Specific Heat Capacity (J/kgK): ');
                disp(' ');

                % Axle Information
                disp('<strong>AXLE INFORMATION:</strong>');
                disp('-----------------');


                % Option to dynamically estimate axle data or enter manually
                axleChoice = MotorThermalAnalysis.getYorNInput('Do you want to dynamically estimate axle data? (Y/N): ');


                if upper(axleChoice) == 'Y'
                    % Logic to dynamically estimate axle data
                    % Example: Estimating based on rotor/stator dimensions
                    axle.OuterDiameter = rotor.InnerDiameter; % Example estimation
                    axle.InnerDiameter = 0; % Typically zero for a solid axle
                    axle.Length = stator.Length; % Example estimation
                    axle.K_Conductivity = 50; % Default or estimated value
                    axle.C_HeatCapacity = 0.466; % Default or estimated value
                else
                    % Manual input for axle data
                    axle.OuterDiameter = input('Enter Axle Outer Diameter (mm): ');
                    axle.InnerDiameter = input('Enter Axle Inner Diameter (mm): '); % Typically zero
                    axle.Length = input('Enter Axle Length (mm): ');
                    axle.K_Conductivity = input('Enter Axle Thermal Conductivity (W/mK): ');
                    axle.C_HeatCapacity = input('Enter Axle Specific Heat Capacity (J/kgK): ');
                end


                % Ask user to choose between default values or custom values for operational parameters
                choiceOperationalParams = MotorThermalAnalysis.getYorNInput('Do you want to use default values for operational parameters? (Y/N): ');

                if upper(choiceOperationalParams) == 'Y'
                    % Use default operational parameters
                    motorOperationalParams = struct('Voltage', 14.7, 'Current', 10, 'Kv', 910, 'PhaseResistance', 0.055178, 'Tmax',defaultTmax);
                else
                    % Prompt for custom operational parameters

                    motorOperationalParams.Voltage = input('Enter operating voltage (V): ');
                    motorOperationalParams.Current = input('Enter operating current (A): ');
                    motorOperationalParams.Kv = input('Enter motor Kv rating (RPM/V): ');
                    motorOperationalParams.PhaseResistance = input('Enter phase resistance (Ohms): ');
                    motorOperationalParams.Tmax = input('Enter maximum operating temperature of motor (C): ');
                end
                % Construct and return the motorComponents and motorOperationalParams structures
                motorComponents = struct('stator', stator, 'rotor', rotor, 'axle', axle);


            end

        end

        function R_thermal = calcThermalResistance(geometry)
            % Extract geometric properties from structure
            outerDiameter = geometry.OuterDiameter;
            innerDiameter = geometry.InnerDiameter;
            length = geometry.Length;

            % Extract Material Properties from structure
            k = geometry.K_Conductivity;

            % Calculate the thermal Resistance for hollow cylinder
            outerRadius = outerDiameter / 2;
            innerRadius = innerDiameter / 2;

            R_thermal = log(outerRadius/innerRadius) / (2 * pi * k * length);

            % Debugging print statement
            %disp(['<strong>R_cond calculated: </strong>', num2str(R_thermal)]);
        end

        function response = getYorNInput(promptMessage)
            validInput = false;
            while ~validInput
                response = input(promptMessage, 's');
                if strcmpi(response, 'Y') || strcmpi(response, 'N')
                    validInput = true;
                else
                    disp('Invalid input. Please enter Y (Yes) or N (No).');
                end
            end
        end

        function axle = estimateAxleDimensions(stator, rotor)
            % This logic can be adjusted based on how the axle fits with the rotor and stator

            % Estimate the outer diameter of the axle (could be slightly less than the inner diameter of the rotor)
            axle.OuterDiameter = rotor.InnerDiameter * 0.98; % Example estimation

            % The length of the axle could be estimated based on the length of the stator or rotor
            axle.Length = stator.Length; % Example estimation

            % Assign default or estimated values for material properties
            axle.K_Conductivity = 50; % Default or estimated value
            axle.C_HeatCapacity = 0.466; % Default or estimated value

            % Inner diameter is typically zero for a solid axle
            axle.InnerDiameter = 0;
        end

        function Q = estimateHeatTransfer(axle, stator, T_axle, T_stator)
            % Calculate the cross-sectional area of the axle
            A_cross_section = pi * (axle.OuterDiameter / 2)^2;

            % Calculate the length of the heat path
            L_heat_path = axle.Length;

            % Calculate the temperature gradient
            delta_T = T_axle - T_stator;

            % Fourier's law of heat conduction
            Q = axle.K_Conductivity * A_cross_section * delta_T / L_heat_path;

            % Debugging print statement
            %disp(['<strong>Q (Heat Transfer via Axle) calculated: </strong>', num2str(Q)]);
        end

        function [motorSpeed, powerLoss] = estimateMotorParameters(V, I, Kv, R_phase)
            % Adjust motor speed based on current, with a maximum speed of 8000 RPM
            if I < 2
                motorSpeed = 0; % Motor speed approaches zero for low current
            else
                % Scaling factor to adjust motor speed
                scalingFactor = 8000 / (Kv * V); % Adjust this based on maximum RPM and motor characteristics

                % Calculate motor speed with a cap at 8000 RPM
                motorSpeed = min(Kv * V * I * scalingFactor, 8000);
            end

            % Power loss calculation
            powerLoss = I^2 * R_phase;
        end


        function h = calcConvectiveHeatTransferCoeff(motorSpeed, Diameter, T_surface, T_ambient)
            % Calculate the convective heat transfer coefficient (h)
            % motorSpeed in RPM, rotorDiameter in mm, V_inf in m/s

            % Air properties at room temperature (20°C)
            k_air = 0.0253; % Thermal conductivity of air (W/mK)
            Pr = 0.725;     % Prandtl number for air
            beta = 1/293;   % Thermal expansion coefficient (1/K) (assumed for air at ~20°C)
            g = 9.81;       % Acceleration due to gravity (m/s^2)
            nu = 1.568e-5;  % Kinematic viscosity of air (m^2/s)

            % Convert RPM to m/s (assuming the tip of the rotor)
            tip_speed = (motorSpeed * pi / 30) * (Diameter / 1000 / 2);

            % Calculate Reynolds number based on V_inf and rotor diameter
            % (Determine air density and dynamic viscosity at average temperature)
            p = 1.12105; % Density of air
            u = 1.1916e-5; % Dynamic viscosity of air
            Re = (tip_speed * (Diameter / 1000) * p) / u;

            % Define the correlation table
            correlation_table = struct(...
                'Re_ranges', {[0.4, 4; 4, 40; 40, 4000; 4000, 40000; 40000, 400000]},...
                'B_values', [0.989, 0.911, 0.683, 0.193, 0.027],...
                'n_values', [0.330, 0.385, 0.366, 0.618, 0.805]...
                );
            isForcedConvectionApplicable = false;
            for i = 1:size(correlation_table.Re_ranges, 1)
                if Re >= correlation_table.Re_ranges(i, 1) && Re < correlation_table.Re_ranges(i, 2)
                    B = correlation_table.B_values(i);
                    n = correlation_table.n_values(i);
                    isForcedConvectionApplicable = true;
                    break;
                end
            end

            % Calculate h based on the type of convection
            if isForcedConvectionApplicable
                % Forced Convection
                Nu = B * Re^n * Pr^(1/3);
            else
                % Natural Convection
                delta_T = T_surface - T_ambient;
                Gr = g * beta * abs(delta_T) * (Diameter / 1000)^3 / nu^2; % Grashof number
                Ra = Gr * Pr; % Rayleigh number

                % Empirical correlation for natural convection (e.g., vertical plate)
                if Ra < 1e9
                    Nu = 0.59 * Ra^(1/4);
                else
                    Nu = 0.10 * Ra^(1/3);
                end
            end

            % Solve for h (using rotor diameter as characteristic length)
            h = Nu * k_air / (Diameter / 1000); % h in W/m^2K
        end

        function R_conv = calcConvectiveResistance(h, surfaceArea)
            % Calculate the equivalent convective resistance (R_conv)
            % h is the convective heat transfer coefficient (W/m^2K)
            % surfaceArea in m^2

            R_conv = 1 / (h * surfaceArea);
            %fprintf("Here is the R_conv: %d",R_conv);
        end

        function I = currentProfile(t, lowCurrentStarts, lowCurrentEnds)
            % Simplified current profile for constant RPM
            % Check if it's a low-current period
            if MotorThermalAnalysis.isLowCurrentPeriod(t, lowCurrentStarts, lowCurrentEnds)
                I = 1; % Low current value
                return;
            end

            % Existing logic for normal operation
            % Base Oscillation
            I = 15 + 6 * sin(2 * pi * t / 50);
        end

        function isLowCurrent = isLowCurrentPeriod(t, lowCurrentStarts, lowCurrentEnds)
            % Check if current time is within any of the low current periods
            isLowCurrent = any(t >= lowCurrentStarts & t <= lowCurrentEnds);
            %disp(['Low Current Period Check at time: ', num2str(t), ', Result: ', num2str(isLowCurrent)]);
        end

        function I = currentProfileForGait(t, total_time, lowCurrentStarts, lowCurrentEnds)
            % Check for low current period
            if MotorThermalAnalysis.isLowCurrentPeriod(t, lowCurrentStarts, lowCurrentEnds)
                I = 1; % Low current value for low current periods
                return;
            end

            % Divide the total time into three phases: walking, trotting, and running
            phase_duration = total_time / 3;

            % Walking Phase (0 to phase_duration)
            if t <= phase_duration
                % Moderate and stable current for walking
                I = 5 + 1 * sin(2 * pi * t / 50); % Example values

                % Trotting Phase (phase_duration to 2 * phase_duration)
            elseif t <= 2 * phase_duration
                % Higher and more fluctuating current for trotting
                I = 8 + 3 * sin(2 * pi * t / 25); % More fluctuation

                % Running Phase (2 * phase_duration to total_time)
            else
                % Highest and most fluctuating current for running
                I = 15 + 8 * sin(2 * pi * t / 15); % Highest peaks
            end
        end
    end

end