%% Stride length on simulated data
% Zahraa Sattar | EEE4022S
%
% Integrate the torso-velocity track between detected step events.
% Step recognition and walking-speed results are taken from the existing workspace.

torsoVelocityMethod = 'walkingspeed';

if strcmp(torsoVelocityMethod, 'walkingspeed')   % swap in the more refined walking-speed pipeline instead
    
    torsoVelocity_m_s = velocityTrack_mps;
    torsoVelocity_m_s = fillmissing(torsoVelocity_m_s, 'linear', 'EndValues', 'nearest');
    torsoVelocityTime_s = final.time_axis;              % matching time axis from that script's output

else   

    [S_torso, freq_axis_torso, time_axis_torso] = Zahraa_Spectrogram(RdrRxSignal, win, noverlap, nfft, fs);

    torsoBand_Hz = [5, 350*fc/5.8e9];   % widened band at 350 Hz, restricts search to torso return bandwidth

    bandIdx = find(freq_axis_torso >= torsoBand_Hz(1) & freq_axis_torso <= torsoBand_Hz(2));    % frequency-axis indices that fall inside the torso band

    torsoFreq_Hz = zeros(1, size(S_torso,2));                   % one torso frequency value per time frame, filled in below
    
    for tIdx = 1:size(S_torso,2)                                % loop over every time frame in the spectrogram
        P_col = abs(S_torso(bandIdx,tIdx)).^2;                  % power in the torso band, this frame only
        [~, relIdx] = max(P_col);
        torsoFreq_Hz(tIdx) = freq_axis_torso(bandIdx(relIdx));
    end

    torsoVelocity_m_s = torsoFreq_Hz*lamda/2;   % v = f*lambda/2 (standard Doppler equation)
    torsoVelocityTime_s = time_axis_torso;      % keeps the matching time axis alongside the velocity track

end

%% Step 3: Ground truth (heel-strike based) - correct reference for stride LENGTH

trueStepTimes_s = bufferTime_s : (Tc/2) : (bufferTime_s + SimulationTime_s);

%% Step 4: Ground-truth stride length, from the simulator's own torso trajectory
[groundTruth_stepLengths_m, groundTruth_strideLengths_m] = ComputeGroundTruthStrideLength(t, Torso_R, trueStepTimes_s);

%% Step 5: Estimated stride length, from integrated Doppler velocity + detected steps
[estimated_stepLengths_m, estimated_strideLengths_m, displacement_m] = EstimateStrideLength(torsoVelocity_m_s, torsoVelocityTime_s, detectedStepTimes_s);

%% Step 6: Compare
fprintf('--- Stride length comparison ---\n');
fprintf('Ground truth: mean %.4f m, std %.4f m (n=%d)\n', mean(groundTruth_strideLengths_m), std(groundTruth_strideLengths_m), numel(groundTruth_strideLengths_m));     % summary stats for the true stride lengths
fprintf('Estimated:    mean %.4f m, std %.4f m (n=%d)\n', mean(estimated_strideLengths_m), std(estimated_strideLengths_m), numel(estimated_strideLengths_m));           % same summary stats for the estimated stride lengths
fprintf('Constant speed check, v0*Tc = %.4f m\n', v0*Tc);

%% Local functions

function [stepLengths_m, strideLengths_m] = ComputeGroundTruthStrideLength(t, Torso_R, stepTimes_s)

    Torso_R_atSteps = interp1(t, Torso_R, stepTimes_s, 'linear');               % reads the simulator's real torso position at each step time
    stepLengths_m = abs(diff(Torso_R_atSteps));                                 % distance moved between each pair of consecutive steps
    strideLengths_m = abs(Torso_R_atSteps(3:end) - Torso_R_atSteps(1:end-2));   % distance moved between every other step (same foot, one full stride)
end

function [stepLengths_m, strideLengths_m, displacement_m] = EstimateStrideLength(velocity_m_s, velocityTime_s, stepTimes_s)

    displacement_m = cumtrapz(velocityTime_s, velocity_m_s);                                % integrates velocity over time to reconstruct a position track
    displacement_atSteps = interp1(velocityTime_s, displacement_m, stepTimes_s, 'linear');  % reads that reconstructed position at each detected step time
    stepLengths_m = abs(diff(displacement_atSteps));                                        % distance moved between each pair of consecutive detected steps
    strideLengths_m = abs(displacement_atSteps(3:end) - displacement_atSteps(1:end-2));     % distance moved between every other detected step
end
