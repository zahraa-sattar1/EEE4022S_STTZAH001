%% Elbow flexion angle
% Zahraa Sattar | EEE4022S
%
% Calculate elbow flexion in degrees over the gait cycle.
% The walking-speed branch sets the control points, with PCHIP interpolation between phases.

function elbowFlexion_deg = computeElbowFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
if gaitCharacteristic == 'a'
    elbow_x1 = 0.05; elbow_x2 = 0.5; elbow_x3 = 0.9;
    elbow_y1 = 6*relativeWalkingVelocity + 3;
    elbow_y2 = 34*relativeWalkingVelocity + 3;
    elbow_y3 = 10*relativeWalkingVelocity + 3;
elseif gaitCharacteristic == 'b'
    elbow_x1 = 0.05;
    elbow_x2 = 0.01*(relativeWalkingVelocity-0.5)/0.8 + 0.5;
    elbow_x3 = 0.9;
    elbow_y1 = 8*(relativeWalkingVelocity-0.5)/0.8 + 6;
    elbow_y2 = 24*(relativeWalkingVelocity-0.5)/0.8 + 20;
    elbow_y3 = 9*(relativeWalkingVelocity-0.5)/0.8 + 8;
else
    elbow_x1 = 0.05;
    elbow_x2 = 0.04*(relativeWalkingVelocity-1.3)/1.7 + 0.51;
    elbow_x3 = -0.1*(relativeWalkingVelocity-1.3)/1.7 + 0.9;
    elbow_y1 = -6*(relativeWalkingVelocity-1.3)/1.7 + 14;
    elbow_y2 = 26*(relativeWalkingVelocity-1.3)/1.7 + 44;
    elbow_y3 = -6*(relativeWalkingVelocity-1.3)/1.7 + 17;

end
elbowPhasePoints = [elbow_x1-1,elbow_x2-1,elbow_x3-1, ...
                    elbow_x1,elbow_x2,elbow_x3, ...
                    elbow_x1+1,elbow_x2+1,elbow_x3+1];
elbowAnglePoints_deg = [elbow_y1,elbow_y2,elbow_y3, ...
                        elbow_y1,elbow_y2,elbow_y3, ...
                        elbow_y1,elbow_y2,elbow_y3];
elbowFlexion_deg = pchip(elbowPhasePoints,elbowAnglePoints_deg,tR);
end
