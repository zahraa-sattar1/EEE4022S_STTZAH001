%% Measured walking demonstration
% Zahraa Sattar | EEE4022S
%
% Run the ten supplied walking recordings and summarise the four parameters.
% Compare the estimates with the retained results from the completed analysis.

function resultsFolder = Run_Measured_Demo()
root = fileparts(mfilename('fullpath'));
oldPath = path;
pathCleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root, 'Extraction Algorithms on Measured Data'));
resultsFolder = Run_Measured_Recordings();
Summarise_Measured_Results(resultsFolder);
T = readtable(fullfile(resultsFolder, 'Results.csv'), 'TextType', 'string');
E = readtable(fullfile(root, 'Participant Recordings', 'Expected_Results.csv'), 'TextType', 'string');
[found, index] = ismember(T.wav_file, E.wav_file);
assert(all(found) && height(T) == height(E), 'Expected-results recording list differs.');
fields = {'estSpeed_mps','estStepCount','estStrideDuration_mean_s','estStrideLength_mean_m'};
check = table(T.wav_file, 'VariableNames', {'wav_file'});
for k = 1:numel(fields)
    field = fields{k};
    check.([field '_difference']) = T.(field)-E.(field)(index);
end
check.withinTolerance = all(abs(check{:, 2:end}) <= 1e-6, 2);
writetable(check, fullfile(resultsFolder, 'Reproduction_Check.csv'));
fprintf('%d/%d recordings match the retained estimates within 1e-6.\n', sum(check.withinTolerance), height(check));
assert(all(check.withinTolerance), 'Demo:ReproductionMismatch', ...
    'One or more estimates differ. Inspect Reproduction_Check.csv.');
end
