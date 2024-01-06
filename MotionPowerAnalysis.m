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
            %fprintf("Here is force per leg: %d\n",forcePerLeg)

            % Calculate the vertical component of the force assuming the force acts at the end of the link
            verticalForceComponent = forcePerLeg * cos(jointAngleRad);

            % Calculate the torque at the joint Torque
            torque = verticalForceComponent * linkLength;
            %fprintf("Here is the torque for the leg %d",torque);
        end

        function [theta1_up, theta2_up, theta1_down, theta2_down] = calculateJointAnglesFor2D(L1, L2, X, Y)
            % Calculate the distance from the hip to the foot
            D = sqrt(X^2 + Y^2);

            % Check if the point is reachable
            if D > L1 + L2
                error('The point is not reachable by the leg');
            end

            % Calculate the angle from the horizontal to the foot
            gamma = atan2(Y, X);

            % Calculate the Knee Angle for both elbow-up and elbow-down configurations
            cos_alpha = (L1^2 + L2^2 - D^2)/(2*L1*L2);
            alpha_up = acos(cos_alpha);
            alpha_down = -alpha_up;

            % Calculate the Hip Angle for both configurations
            cos_beta = (D^2 + L1^2 - L2^2)/(2*D*L1);
            beta_up = acos(cos_beta);
            beta_down = -beta_up;

            % Calculate the Hip Joint angles for both configurations
            theta2_up = gamma - beta_up;
            theta2_down = gamma - beta_down;

            % The Knee angles are complementary to alpha
            theta1_up = pi - alpha_up;
            theta1_down = pi - alpha_down;

            % Convert angles to degrees
            theta1_up = rad2deg(theta1_up);
            theta2_up = rad2deg(theta2_up);
            theta1_down = rad2deg(theta1_down);
            theta2_down = rad2deg(theta2_down);
        end

        %function [amax, jmax] = EstimateMotorCapabilities(torque, insertia)
        function [rotorInertia, axleInertia] = calculateMotorInertia(outerDiameterRotor, innerDiameterRotor, lengthRotor, densityRotor, outerDiameterAxle, lengthAxle, densityAxle)
            % Rotor First
            outerRadiusRotor = outerDiameterRotor / 2;
            innerRadiusRotor = innerDiameterRotor / 2;

            volumeRotor = pi * (outerRadiusRotor^2 - innerRadiusRotor^2) * lengthRotor;
            massRotor = volumeRotor * densityRotor;

            rotorInertia = 0.5 * massRotor * outerRadiusRotor^2;

            % Axle
            radiusAxle = outerDiameterAxle / 2;
            volumeAxle = pi * radiusAxle^2 * lengthAxle;
            massAxle = volumeAxle * densityAxle;
            axleInertia = 0.5 * massAxle * radiusAxle^2;
        end

        function totalLoadInertia = calculateLoadInertia(massUpperLink, lengthUpperLink, massLowerLink, lengthLowerLink, angleLowerLink)
            % Calculate the inertia of the upper link about its center of mass
            % Assuming each link can be approximated as a rod pivoted at one end
            inertiaUpperLink = (1/3) * massUpperLink * lengthUpperLink^2;
            disp(['distal (upper) link: ', num2str(inertiaUpperLink)]);

            % Calculate the vertical distance from the hip to the COM of the lower link
            % Assuming the lower link rotates about the end of the upper link
            distanceLowerLink = (lengthUpperLink + lengthLowerLink/2) * sind(angleLowerLink);
            disp(['distance to hip joint: ', num2str(distanceLowerLink)]);

            % Calculate the inertia of the lower link about its center of mass
            inertiaLowerLinkCenter = (1/12) * massLowerLink * lengthLowerLink^2;
            disp(['proximal (lower) link: ', num2str(inertiaLowerLinkCenter)]);

            % Use the parallel axis theorem to find the inertia of the lower link about the hip
            inertiaLowerLink = inertiaLowerLinkCenter + massLowerLink * distanceLowerLink^2;
            disp(['proximal (lower) link total: ', num2str(inertiaLowerLink)]);

            % Sum the inertias to get the total load inertia
            totalLoadInertia = (inertiaUpperLink + inertiaLowerLink);
            disp(['total proximal link: ', num2str(totalLoadInertia)]);
        end

        function [q, v, acc, j, t] = generateSCurveTrajectory(amax, vmax, qf, jmaxRatio)

            %% 7-Segement S Curve Trajectory Planning
            jmax = jmaxRatio * amax; % Adjust jmax based on amax since I have no good way to estimate a maximum desired impulse since it based on teh construction of the acuator and the leg themselves
            %jmax=12; % Maximum Jerk the system can tolerate (comes from the physical system)
            %amax=3; % Maximum acceleration the actuator can produce
            %vmax=2; % Maxium velocity the actuator can produce
            T = amax/jmax + vmax/amax + qf/vmax; % Total time required

            % Using the EQs garned from the 4 additional conditions we have left over
            t1 = amax/jmax;
            t3 = vmax/amax + t1;
            t2 = t3 - t1;
            % Since we have symmetry across the midline we can state the following
            t4 = T - t3;
            t5 = T - t2;
            t6 = T - t1;


            %% Next we need to solve for the unknown coefficients
            % Angle q is assumed to be of cubic form and has 4 coef
            % Velocity v therefore is of quadratic form and has 3 coef
            % Acceleration a therefore is of linear form and has 2 coef
            % Jerk is therefore a constant value and has 1 coef

            % using symmetry, continutity conditions, and clever inital and final
            % conditions we can reduce the 28 unknown coef down to 14
            a1 = jmax;
            a2 = 0;
            a3 = -jmax;
            a4 = a2;
            a5 = a3;
            a6 = a2;
            a7 = a1;

            b1 = 0;
            b2 = amax;
            b3 = amax + jmax*t2;
            b4 = 0;
            b5 = jmax*t4;
            b6 = -amax;
            b7 = -amax-jmax*t6;

            c1 = 0;
            c2 = ((a1*t1^2)/2+b1*t1+c1)-((a2*t1^2)/2+b2*t1);
            c3 = ((a2*t2^2)/2+b2*t2+c2)-((a3*t2^2)/2+b3*t2);
            c4 = ((a3*t3^2)/2+b3*t3+c3)-((a4*t3^2)/2+b4*t3);
            c5 = ((a4*t4^2)/2+b4*t4+c4)-((a5*t4^2)/2+b5*t4);
            c6 = ((a5*t5^2)/2+b5*t5+c5)-((a6*t5^2)/2+b6*t5);
            c7 = ((a6*t6^2)/2+b6*t6+c6)-((a7*t6^2)/2+b7*t6);

            d1 = 0;
            d2 = ((a1*t1^3)/6+(b1*t1^2)/2+c1*t1+d1)-((a2*t1^3)/6+(b2*t1^2)/2+c2*t1);
            d3 = ((a2*t2^3)/6+(b2*t2^2)/2+c2*t2+d2)-((a3*t2^3)/6+(b3*t2^2)/2+c3*t2);
            d4 = ((a3*t3^3)/6+(b3*t3^2)/2+c3*t3+d3)-((a4*t3^3)/6+(b4*t3^2)/2+c4*t3);
            d5 = ((a4*t4^3)/6+(b4*t4^2)/2+c4*t4+d4)-((a5*t4^3)/6+(b5*t4^2)/2+c5*t4);
            d6 = ((a5*t5^3)/6+(b5*t5^2)/2+c5*t5+d5)-((a6*t5^3)/6+(b6*t5^2)/2+c6*t5);
            d7 = ((a6*t6^3)/6+(b6*t6^2)/2+c6*t6+d6)-((a7*t6^3)/6+(b7*t6^2)/2+c7*t6);

            % Establish time vector
            N=100;
            t=linspace(0,T,N);

            ind1 = find((0<=t) & (t<=t1));
            ind2 = find((t1<=t) & (t<=t2));
            ind3 = find((t2<=t) & (t<=t3));
            ind4 = find((t3<=t) & (t<=t4));
            ind5 = find((t4<=t) & (t<=t5));
            ind6 = find((t5<=t) & (t<=t6));
            ind7 = find((t6<=t) & (t<=T));

            % Establish Angle Equations
            q1 = a1*t.^3/6+b1*t.^2/2+c1*t+d1;
            q2 = a2*t.^3/6+b2*t.^2/2+c2*t+d2;
            q3 = a3*t.^3/6+b3*t.^2/2+c3*t+d3;
            q4 = a4*t.^3/6+b4*t.^2/2+c4*t+d4;
            q5 = a5*t.^3/6+b5*t.^2/2+c5*t+d5;
            q6 = a6*t.^3/6+b6*t.^2/2+c6*t+d6;
            q7 = a7*t.^3/6+b7*t.^2/2+c7*t+d7;
            q=[q1(ind1) q2(ind2) q3(ind3) q4(ind4) q5(ind5) q6(ind6) q7(ind7)];

            % Establish Velocity Equations
            v1 = a1*t.^2/2+b1*t+c1;
            v2 = a2*t.^2/2+b2*t+c2;
            v3 = a3*t.^2/2+b3*t+c3;
            v4 = a4*t.^2/2+b4*t+c4;
            v5 = a5*t.^2/2+b5*t+c5;
            v6 = a6*t.^2/2+b6*t+c6;
            v7 = a7*t.^2/2+b7*t+c7;
            v=[v1(ind1) v2(ind2) v3(ind3) v4(ind4) v5(ind5) v6(ind6) v7(ind7)];

            % Establish Acceleration Equations
            acc1=a1*t+b1;
            acc2=a2*t+b2;
            acc3=a3*t+b3;
            acc4=a4*t+b4;
            acc5=a5*t+b5;
            acc6=a6*t+b6;
            acc7=a7*t+b7;
            acc=[acc1(ind1) acc2(ind2) acc3(ind3) acc4(ind4) acc5(ind5) acc6(ind6) acc7(ind7)];

            % Establish Jerk Equations
            j1=a1*ones(1,length(t));
            j2=a2*ones(1,length(t));
            j3=a3*ones(1,length(t));
            j4=a4*ones(1,length(t));
            j5=a5*ones(1,length(t));
            j6=a6*ones(1,length(t));
            j7=a7*ones(1,length(t));
            j=[j1(ind1) j2(ind2) j3(ind3) j4(ind4) j5(ind5) j6(ind6) j7(ind7)];

            %% Plotting the Results
            figure(4);
            % Plot Angle
            subplot(4,1,1);
            plot(t,q,'b','LineWidth',2);
            xlim([0,T]);
            xlabel('t (sec)');
            ylabel('q(t) (rad)');
            title('S curve (7-Segment)');

            % Plot Velocity
            subplot(4,1,2);
            plot(t,v,'g','LineWidth',2);
            xlim([0,T]);
            xlabel('t (sec)');
            ylabel('dq/dt (t) (rad/s)');

            % Plot Acceleration
            subplot(4,1,3);
            plot(t,acc,'r','LineWidth',2);
            xlim([0,T]);
            xlabel('t (sec)');
            ylabel('d^[2]/dt^[2] (t) (rad/s^[2])');

            % Plot Jerk
            subplot(4,1,4);
            plot(t,j,'k','LineWidth',2);
            xlim([0,T]);
            xlabel('t (sec)');
            ylabel('d^[3]/dt^[3] (t) (rad/s^[3])');

        end

        function [z_position, x_position] = calculateLegTrajectory(t, Rh, H, Sx, Ts, Ty, T)
            % Ensure that t is within the range [0, T]
            if t < 0 || t > T
                error('Time t must be within the gait period [0, T].');
            end

            % Calculate alpha and tau
            alpha = t - Ts/2;
            tau = alpha / Ty;

            % Calculate the vertical position z(t)
            if t <= Ts/2
                z_position = -Rh;
            elseif alpha <= Ty/2 && t <= Ty + Ts/2
                z_position = 2*H*(tau - (1/(4*pi))*sin(4*pi*tau)) - Rh;
            elseif alpha > Ty/2 && t <= Ty + Ts/2
                z_position = -2*H*(tau - (1/(4*pi))*sin(4*pi*tau)) - Rh;
            else
                z_position = -Rh;
            end

            % Calculate the horizontal position x(t)
            if t <= Ts/2
                x_position = -(Sx/Ts)*t;
            elseif t <= Ty + Ts/2
                x_position = Sx*(tau - (sin(2*pi*tau))/(2*pi)) - Sx/2;
            else
                x_position = 0.5*Sx*tau + 0.5*Sx*(1 - ((Ty + 0.5*Ts - tau)/(Ty + 0.5*Ts - T)));
            end
        end
    end
end