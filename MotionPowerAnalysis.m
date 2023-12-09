classdef MotionPowerAnalysis
    methods (Static)
        % Function to plot a solid cylinder
        function plotSolidCylinder(outerDiameter, innerDiameter, length, color, zOffset)
            numPoints = 100; % Define the number of points on the circumference
            [X, Y, Z] = cylinder([outerDiameter/2, outerDiameter/2], numPoints);
            surf(X, Y, Z*length + zOffset, 'FaceColor', color); % Outer cylinder

            if innerDiameter > 0
                [Xi, Yi, Zi] = cylinder([innerDiameter/2, innerDiameter/2], numPoints);
                surf(Xi, Yi, Zi*length + zOffset, 'FaceColor', 'white'); % Inner cylinder (hollow part)
            end

            % Adding top and bottom caps correctly
            theta = linspace(0, 2*pi, numPoints);
            [XOuterCap, YOuterCap] = pol2cart(theta, outerDiameter/2);
            [XInnerCap, YInnerCap] = pol2cart(theta, innerDiameter/2);

            nCaps = size(XOuterCap, 2) - 1; % Number of segments for the cap

            % Top cap
            for i = 1:nCaps
                fill3([XOuterCap(i:i+1), XInnerCap(i+1:-1:i)], ...
                    [YOuterCap(i:i+1), YInnerCap(i+1:-1:i)], ...
                    [length + zOffset, length + zOffset, length + zOffset, length + zOffset], color);
            end

            % Bottom cap
            for i = 1:nCaps
                fill3([XOuterCap(i:i+1), XInnerCap(i+1:-1:i)], ...
                    [YOuterCap(i:i+1), YInnerCap(i+1:-1:i)], ...
                    [zOffset, zOffset, zOffset, zOffset], color);
            end
        end
    end
end