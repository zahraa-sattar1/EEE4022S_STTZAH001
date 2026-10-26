%% Ankle flexion angle
% Zahraa Sattar | EEE4022S
%
% Calculate ankle flexion in degrees over the gait cycle.
% The walking-speed branch sets the control points, with PCHIP interpolation between phases.

function ankleFlexion_deg = computeAnkleFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,relativeSupportDuration,tR)
if gaitCharacteristic == 'a'
    ankle_x1 = 0;
    ankle_x2 = 0.08;
    ankle_x3 = 0.5;
    ankle_x4 = relativeSupportDuration;
    ankle_x5 = 0.85;
    ankle_y1 = -3;
    ankle_y2 = -30*relativeWalkingVelocity - 3;
    ankle_y3 = 22*relativeWalkingVelocity - 3;
    ankle_y4 = -34*relativeWalkingVelocity - 3;
    ankle_y5 = -3;
elseif gaitCharacteristic == 'b'
    ankle_x1 = 0;
    ankle_x2 = 0.08;
    ankle_x3 = 0.5;
    ankle_x4 = relativeSupportDuration;
    ankle_x5 = 0.85;
    ankle_y1 = -3;
    ankle_y2 = -18;
    ankle_y3 = 8;
    ankle_y4 = -20;
    ankle_y5 = -3;
else
    ankle_x1 = 0;
    ankle_x2 = 0.08;
    ankle_x3 = -0.1*(relativeWalkingVelocity-1.3)/1.7 + 0.5;
    ankle_x4 = relativeSupportDuration;
    ankle_x5 = 0.85;
    ankle_y1 = 5*(relativeWalkingVelocity-1.3)/1.7 - 3;
    ankle_y2 = 4*(relativeWalkingVelocity-1.3)/1.7 - 18;
    ankle_y3 = -3*(relativeWalkingVelocity-1.3)/1.7 + 8;
    ankle_y4 = -8*(relativeWalkingVelocity-1.3)/1.7 - 20;
    ankle_y5 = 5*(relativeWalkingVelocity-1.3)/1.7 - 3;
end
anklePhasePoints = [ankle_x1-1,ankle_x2-1,ankle_x3-1,ankle_x4-1,ankle_x5-1, ...
                    ankle_x1,ankle_x2,ankle_x3,ankle_x4,ankle_x5, ...
                    ankle_x1+1,ankle_x2+1,ankle_x3+1,ankle_x4+1,ankle_x5+1];
ankleAnglePoints_deg = [ankle_y1,ankle_y2,ankle_y3,ankle_y4,ankle_y5, ...
                        ankle_y1,ankle_y2,ankle_y3,ankle_y4,ankle_y5, ...
                        ankle_y1,ankle_y2,ankle_y3,ankle_y4,ankle_y5];
ankleFlexion_deg = pchip(anklePhasePoints,ankleAnglePoints_deg,tR);
end
