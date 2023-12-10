classdef MotionPowerAnalysis
    methods (Static)
        % Function to plot a solid cylinder
        function plotSolidCylinder(ax, outerDiameter, innerDiameter, length, color, zOffset)
            %disp(['Plotting cylinder with outer diameter: ', num2str(outerDiameter)]);
            %disp(['Inner diameter: ', num2str(innerDiameter)]);
            %disp(['Length: ', num2str(length)]);
            %disp(['Color: ', color]);
            %disp(['zOffset: ', num2str(zOffset)]);

            numPoints = 100; % Define the number of points on the circumference
            [X, Y, Z] = cylinder([outerDiameter/2, outerDiameter/2], numPoints);

            % Debug: Size of the matrices X, Y, Z
            %disp(['Size of X: ', mat2str(size(X))]);
            %disp(['Size of Y: ', mat2str(size(Y))]);
            %disp(['Size of Z: ', mat2str(size(Z))]);

            surf(ax, X, Y, Z*length + zOffset, 'FaceColor', color);

            % Create the inner cylinder (if inner diameter is specified)
            if innerDiameter > 0
                [Xi, Yi, Zi] = cylinder([innerDiameter/2, innerDiameter/2], numPoints);
                surf(ax, Xi, Yi, Zi*length + zOffset, 'FaceColor', 'white');
            end

            % Add top and bottom caps
            theta = linspace(0, 2*pi, numPoints);
            [XOuterCap, YOuterCap] = pol2cart(theta, outerDiameter/2);
            [XInnerCap, YInnerCap] = pol2cart(theta, innerDiameter/2);

            nCaps = size(XOuterCap, 2) - 1; % Number of segments for the cap
            %disp(['Number of cap segments: ', num2str(nCaps)]);

            % Top cap
            for i = 1:nCaps
                fill3(ax, [XOuterCap(i:i+1), XInnerCap(i+1:-1:i)], ...
                      [YOuterCap(i:i+1), YInnerCap(i+1:-1:i)], ...
                      [length + zOffset, length + zOffset, length + zOffset, length + zOffset], color);
            end

            % Bottom cap
            for i = 1:nCaps
                fill3(ax, [XOuterCap(i:i+1), XInnerCap(i+1:-1:i)], ...
                      [YOuterCap(i:i+1), YInnerCap(i+1:-1:i)], ...
                      [zOffset, zOffset, zOffset, zOffset], color);
            end
        end
    end
end