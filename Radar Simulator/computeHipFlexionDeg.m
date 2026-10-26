%% Hip flexion angle
% Zahraa Sattar | EEE4022S
%
% Calculate hip flexion in degrees over the gait cycle.
% The walking-speed branch sets the control points, with PCHIP interpolation between phases.

function hipFlexion_deg = computeHipFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
if gaitCharacteristic == 'a'
    hip_x1 = -0.1; hip_x2 = 0.5; hip_x3 = 0.9;
    hip_y1 = 50*relativeWalkingVelocity;
    hip_y2 = -30*relativeWalkingVelocity;
    hip_y3 = 50*relativeWalkingVelocity;
elseif gaitCharacteristic == 'b'
    hip_x1 = -0.1; hip_x2 = 0.5; hip_x3 = 0.9;
    hip_y1 = 25; hip_y2 = -15; hip_y3 = 25;
else
    hip_x1 = 0.2*(relativeWalkingVelocity-1.3)/1.7 - 0.1;
    hip_x2 = 0.5;
    hip_x3 = 0.9;
    hip_y1 = 5*(relativeWalkingVelocity-1.3)/1.7 + 25;
    hip_y2 = -15;
    hip_y3 = 6*(relativeWalkingVelocity-1.3)/1.7 + 25;
end
if abs((hip_x1 + 1) - hip_x3) < 1e-12
    hipPhasePoints = [hip_x1-1,hip_x2-1,hip_x1,hip_x2,hip_x3,hip_x2+1,hip_x3+1];
    hipAnglePoints_deg = [hip_y1,hip_y2,hip_y1,hip_y2,hip_y3,hip_y2,hip_y3];
else
    hipPhasePoints = [hip_x1-1,hip_x2-1,hip_x3-1,hip_x1,hip_x2,hip_x3, ...
                      hip_x1+1,hip_x2+1,hip_x3+1];
    hipAnglePoints_deg = [hip_y1,hip_y2,hip_y3,hip_y1,hip_y2,hip_y3, ...
                          hip_y1,hip_y2,hip_y3];
end
hipFlexion_deg = pchip(hipPhasePoints,hipAnglePoints_deg,tR);
end
