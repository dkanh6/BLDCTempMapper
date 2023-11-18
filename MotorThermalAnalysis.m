
classdef MotorThermalAnalysis
    methods (Static)

        function [cp_composite, k_composite] = calcCompositeProperties(copperFraction)
            % Define the specific heat capacities (J/g°C)
            cp_copper = 0.385; % Copper
            cp_steel = 0.466;  % Steel

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

        function motorComponents = getUserInput()
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
                rotor = struct('OuterDiameter', 35, 'InnerDiameter', 30, 'Length', 36, 'K_Conductivity', 167, 'C_HeatCapacity', 0.89);
                axle = struct('OuterDiameter', 5, 'InnerDiameter', 0, 'Length', 55, 'K_Conductivity', 50, 'C_HeatCapacity', 0.466);
            else
                % Prompt for custom values
                disp('==============================');
                disp('CUSTOM MOTOR COMPONENT INPUT');
                disp('==============================');

                % Stator Information
                disp('STATOR INFORMATION:');
                disp('-------------------');
                stator.OuterDiameter = input('Enter Stator Outer Diameter (mm): ');
                stator.InnerDiameter = input('Enter Stator Inner Diameter (mm): ');
                stator.Length = input('Enter Stator Length (mm): ');
                stator.K_Conductivity = input('Enter Stator Thermal Conductivity (W/mK): ');
                stator.C_HeatCapacity = input('Enter Stator Specific Heat Capacity (J/kgK): ');
                disp(' ');

                % Rotor Information
                disp('ROTOR INFORMATION:');
                disp('------------------');
                rotor.OuterDiameter = input('Enter Rotor Outer Diameter (mm): ');
                rotor.InnerDiameter = input('Enter Rotor Inner Diameter (mm): ');
                rotor.Length = input('Enter Rotor Length (mm): ');
                rotor.K_Conductivity = input('Enter Rotor Thermal Conductivity (W/mK): ');
                rotor.C_HeatCapacity = input('Enter Rotor Specific Heat Capacity (J/kgK): ');
                disp(' ');

                % Axle Information
                disp('AXLE INFORMATION:');
                disp('-----------------');

                % Option to dynamically estimate axle data or enter manually
                axleChoice = MotorThermalAnalysis.getYorNInput('Do you want to dynamically estimate axle data? (Y/N): ');


                if upper(axleChoice) == 'Y'
                    % Logic to dynamically estimate axle data
                    % Example: Estimating based on rotor/stator dimensions
                    axle.OuterDiameter = rotor.InnerDiameter; % Example estimation
                    axle.InnerDiameter = 0; % Typically zero for a solid axle
                    axle.Length = stator.Length; % Example estimation
                    axle.K_Conductivity = 60; % Default or estimated value
                    axle.C_HeatCapacity = 920; % Default or estimated value
                else
                    % Manual input for axle data
                    axle.OuterDiameter = input('Enter Axle Outer Diameter (mm): ');
                    axle.InnerDiameter = input('Enter Axle Inner Diameter (mm): '); % Typically zero
                    axle.Length = input('Enter Axle Length (mm): ');
                    axle.K_Conductivity = input('Enter Axle Thermal Conductivity (W/mK): ');
                    axle.C_HeatCapacity = input('Enter Axle Specific Heat Capacity (J/kgK): ');
                end

                % Construct and return the motorComponents structure
                motorComponents = struct('stator', stator, 'rotor', rotor, 'axle', axle);
            end

        end

        function R_thermal = calcThermalResistance(geometry)
            % Extract geometric properties from structure
            outerDiameter = geometry.OuterDiameter;
            innerDiameter = geometry.InnerDiameter;
            length  = geometry.Length;

            % Extract Material Properties from structure
            k = geometry.K_Conductivity;
            cp = geometry.C_HeatCapacity;

            % Calculate the thermal Resistance for hollow cylinder
            outerRadius = outerDiameter /2;
            innerRadius = innerDiameter /2;

            R_thermal = log(outerRadius/innerRadius) / (2 * pi * k * length);
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
            % Axle and stator structures are expected to contain material properties and dimensions
            % T_axle and T_stator are the temperatures of the axle and stator, respectively

            % Calculate the cross-sectional area of the axle
            A_cross_section = pi * (axle.OuterDiameter / 2)^2;

            % Calculate the length of the heat path (assuming it's the length of the axle)
            L_heat_path = axle.Length;

            % Calculate the temperature gradient
            delta_T = T_axle - T_stator;

            % Fourier's law of heat conduction
            Q = axle.K_Conductivity * A_cross_section * delta_T / L_heat_path;
        end

        function [motorSpeed, powerLoss] = estimateMotorParameters(V, I, Kv, R_phase)
            % Estimate motor speed (RPM) assuming no load
            motorSpeed = Kv * V;  % Ideal no-load speed in RPM

            % Calculate power loss due to resistance in windings
            powerLoss = I^2 * R_phase;  % Power loss in Watts
        end

    end

end