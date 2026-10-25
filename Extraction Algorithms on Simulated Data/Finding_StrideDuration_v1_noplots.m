% Calculate stride duration from the detected step intervals.
% Run step recognition first.

groundTruth_strideDurations_s = Tc*ones(size(strideDurations_s));   % same length as the estimated array, so the two compare one-to-one

%% Step 3: Estimated - already computed by Finding_StepRecognition_v4_noplots's Step 10
estimated_strideDurations_s = strideDurations_s;   

%% Step 4: Compare
fprintf('--- Stride duration comparison ---\n');
fprintf('Ground truth (Tc): %.4f s (n=%d)\n', Tc, numel(groundTruth_strideDurations_s));
fprintf('Estimated:         mean %.4f s, std %.4f s (n=%d)\n', mean(estimated_strideDurations_s,'omitnan'), std(estimated_strideDurations_s,'omitnan'), numel(estimated_strideDurations_s));   % summary stats for the estimated stride durations

