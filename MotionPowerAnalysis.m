classdef MotionPowerAnalysis
    methods (Static)
        % Function to plot a solid cylinder
        function plotSolidCylinder(ax, outerDiameter, innerDiameter, length, color, zOffset)

            numPoints = 100; % Define the number of points on the circumference
            [X, Y, Z] = cylinder([outerDiameter/2, outerDiameter/2], numPoints);

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

        function plotQuadruped(bodyLength, bodyWidth, bodyHeight, upperLegLength, lowerLegLength, shoulderAngles, kneeAngles)
            % Verify that shoulderAngles and kneeAngles are 1x4 arrays
            if length(shoulderAngles) ~= 4 || length(kneeAngles) ~= 4
                error('Shoulder and knee angles must be 1x4 arrays.');
            end

            % Set the axes limits
            maxDimension = max([bodyLength, bodyWidth, bodyHeight, upperLegLength, lowerLegLength]);
            axisLimit = maxDimension * 1.0; % Adjust this factor as needed for the desired scale

            % Define the body as a 3D box
            bodyVertices = [-bodyLength/2, -bodyWidth/2, 0;
                bodyLength/2, -bodyWidth/2, 0;
                bodyLength/2,  bodyWidth/2, 0;
                -bodyLength/2,  bodyWidth/2, 0;
                -bodyLength/2, -bodyWidth/2, bodyHeight;
                bodyLength/2, -bodyWidth/2, bodyHeight;
                bodyLength/2,  bodyWidth/2, bodyHeight;
                -bodyLength/2,  bodyWidth/2, bodyHeight];

            bodyEdges = [1 2; 2 3; 3 4; 4 1; % Lower body square
                5 6; 6 7; 7 8; 8 5; % Upper body square
                1 5; 2 6; 3 7; 4 8]; % Body edges

            % Create figure
            figure;
            hold on;
            grid on;
            axis equal;
            view(25, 30);


            xlim([-axisLimit, axisLimit]);
            ylim([-axisLimit, axisLimit]);
            zlim([-(2*upperLegLength), axisLimit]); % Assuming z starts from 0 (ground level)

            % Define the radius of the spheres representing the joints
            jointRadius = 0.025;

            % Plot the body
            for i = 1:size(bodyEdges, 1)
                plot3(bodyVertices(bodyEdges(i,:), 1), ...
                    bodyVertices(bodyEdges(i,:), 2), ...
                    bodyVertices(bodyEdges(i,:), 3), 'b-','LineWidth',4);
            end

            % Define leg attachment points on the body
            legAttachPoints = [-bodyLength/2, -bodyWidth/2, 0;
                bodyLength/2, -bodyWidth/2, 0;
                bodyLength/2,  bodyWidth/2, 0;
                -bodyLength/2,  bodyWidth/2, 0];

            % Plot the legs
            for i = 1:4
                % Convert angles from degrees to radians for MATLAB functions
                shoulderAngleRad = deg2rad(shoulderAngles(i));
                kneeAngleRad = deg2rad(kneeAngles(i));

                % Upper leg segment (shoulder joint)
                upperLegStart = legAttachPoints(i, :);
                upperLegEnd = upperLegStart + ...
                    [-upperLegLength * sin(shoulderAngleRad), 0, -upperLegLength * cos(shoulderAngleRad)];
                plot3([upperLegStart(1), upperLegEnd(1)], ...
                    [upperLegStart(2), upperLegEnd(2)], ...
                    [upperLegStart(3), upperLegEnd(3)], 'k-','LineWidth',4);

                % Lower leg segment (knee joint)
                lowerLegStart = upperLegEnd;
                lowerLegEnd = lowerLegStart + ...
                    [-lowerLegLength * sin(kneeAngleRad), 0, -lowerLegLength * cos(kneeAngleRad)];
                plot3([lowerLegStart(1), lowerLegEnd(1)], ...
                    [lowerLegStart(2), lowerLegEnd(2)], ...
                    [lowerLegStart(3), lowerLegEnd(3)], 'k-','LineWidth',4);

                % Plot red sphere at the shoulder joint
                [Xs, Ys, Zs] = sphere; % Generate sphere data
                Xs = Xs * jointRadius + upperLegStart(1);
                Ys = Ys * jointRadius + upperLegStart(2);
                Zs = Zs * jointRadius + upperLegStart(3);
                surf(Xs, Ys, Zs, 'EdgeColor', 'none', 'FaceColor', 'r');

                % Plot red sphere at the knee joint
                [Xk, Yk, Zk] = sphere; % Generate sphere data
                Xk = Xk * jointRadius + lowerLegStart(1);
                Yk = Yk * jointRadius + lowerLegStart(2);
                Zk = Zk * jointRadius + lowerLegStart(3);
                surf(Xk, Yk, Zk, 'EdgeColor', 'none', 'FaceColor', 'r');
            end

            % Set labels
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('3D Wireframe of a Quadruped Robot');

            hold off;
        end

        function torque = calculateJointTorque(mass, numLegsOnGround, linkLength, jointAngle, gravity)
            % mass: Total mass of the quadruped (kg)
            % numLegsOnGround: Number of legs on the ground supporting the weight
            % linkLength: Length of the leg link (m)
            % jointAngle: Angle of the joint (degrees)
            % gravity: Acceleration due to gravity (m/s^2), typically 9.81 m/s^2

            % Convert joint angle from degrees to radians
            jointAngleRad = deg2rad(jointAngle);

            % Calculate the force experienced by each leg
            forcePerLeg = (mass * gravity) / numLegsOnGround;
            fprintf("Here is force per leg: %d\n",forcePerLeg)

            % Calculate the vertical component of the force assuming the force acts at the end of the link
            verticalForceComponent = forcePerLeg * cos(jointAngleRad);

            % Calculate the torque at the joint Torque
            torque = verticalForceComponent * linkLength;
        end 
    end
end