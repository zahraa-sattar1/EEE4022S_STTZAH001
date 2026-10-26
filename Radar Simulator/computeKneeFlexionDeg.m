%% Knee flexion angle
% Zahraa Sattar | EEE4022S
%
% Calculate knee flexion in degrees over the gait cycle.
% The walking-speed branch sets the control points, with PCHIP interpolation between phases.

function kneeFlexion_deg = computeKneeFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
if gaitCharacteristic == 'a'
    knee_x1 = 0.17; knee_x2 = 0.4; knee_x3 = 0.75; knee_x4 = 1;
    knee_y1 = 3;
    knee_y2 = 3;
    knee_y3 = 140*relativeWalkingVelocity;
    knee_y4 = 3;
elseif gaitCharacteristic == 'b'
    knee_x1 = 0.17; knee_x2 = 0.4; knee_x3 = 0.75; knee_x4 = 1;
    knee_y1 = 3; knee_y2 = 3; knee_y3 = 70; knee_y4 = 3;
else
    knee_x1 = -0.05*(relativeWalkingVelocity-1.3)/1.7 + 0.17;
    knee_x2 = 0.4;
    knee_x3 = -0.05*(relativeWalkingVelocity-1.3)/1.7 + 0.75;
    knee_x4 = -0.03*(relativeWalkingVelocity-1.3)/1.7 + 1;
    knee_y1 = 22*(relativeWalkingVelocity-1.3)/1.7 + 3;
    knee_y2 = 3;
    knee_y3 = -5*(relativeWalkingVelocity-1.3)/1.7 + 70;
    knee_y4 = 3*(relativeWalkingVelocity-1.3)/1.7 + 3;
end
kneePhasePoints = [knee_x1-1,knee_x2-1,knee_x3-1,knee_x4-1, ...
                   knee_x1,knee_x2,knee_x3,knee_x4, ...
                   knee_x1+1,knee_x2+1,knee_x3+1,knee_x4+1];
kneeAnglePoints_deg = [knee_y1,knee_y2,knee_y3,knee_y4, ...
                       knee_y1,knee_y2,knee_y3,knee_y4, ...
                       knee_y1,knee_y2,knee_y3,knee_y4];
kneeFlexion_deg = pchip(kneePhasePoints,kneeAnglePoints_deg,tR);
end
