
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
            choice = input('Do you want to use default values? (Y/N): ', 's');

            if upper(choice) == 'Y'
                % Use default values
                % Example default values
                % cp copper = 0.385 J/gC
                % cp Aluminum = 0.89 J/gC
                % cp Steel = 0.466 J/gC
                [cp_stator, k_stator] = MotorThermalAnalysis.calcCompositeProperties(0.5); % 50% copper
                stator = struct('OuterDiameter', 27.5, 'InnerDiameter', 5, 'Length', 26, 'K_Conductivity', k_stator, 'C_HeatCapacity', cp_stator); % Assuming 50:50 ratio of copper to steel to calculate the thermal capacitance and conductivity
                rotor = struct('OuterDiameter', 35, 'InnerDiameter', 30, 'Length', 36, 'K_Conductivity', 45, 'C_HeatCapacity', 850);
                axle = struct('OuterDiameter', 20, 'InnerDiameter', 0, 'Length', 60, 'K_Conductivity', 167, 'C_HeatCapacity', 0.89);
            else
                % Prompt for custom values
                disp('Please enter the stator information:');
                stator.OuterDiameter = input('Enter stator outer diameter (mm): ');
                stator.InnerDiameter = input('Enter stator inner diameter (mm): ');
                stator.Length = input('Enter stator length (mm): ');
                stator.K_Conductivity = input('Enter stator thermal conductivity (W/mK): ');
                stator.C_HeatCapacity = input('Enter stator specific heat capacity (J/kgK): ');

                disp('Please enter the rotor information:');
                rotor.OuterDiameter = input('Enter rotor outer diameter (mm): ');
                rotor.InnerDiameter = input('Enter rotor inner diameter (mm): ');
                rotor.Length = input('Enter rotor length (mm): ');
                rotor.K_Conductivity = input('Enter rotor thermal conductivity (W/mK): ');
                rotor.C_HeatCapacity = input('Enter rotor specific heat capacity (J/kgK): ');

                disp('Please enter the axle information:');
                axle.OuterDiameter = input('Enter axle outer diameter (mm): ');
                axle.InnerDiameter = input('Enter axle inner diameter (mm): ');  % Typically zero
                axle.Length = input('Enter axle length (mm): ');
                axle.K_Conductivity = input('Enter axle thermal conductivity (W/mK): ');
                axle.C_HeatCapacity = input('Enter axle specific heat capacity (J/kgK): ');
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



    end

end