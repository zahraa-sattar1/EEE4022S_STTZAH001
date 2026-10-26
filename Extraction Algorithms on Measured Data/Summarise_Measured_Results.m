%% Summarise measured results
% Zahraa Sattar | EEE4022S
%
% Compare the measured-demo estimates with their recording references.
% Save error statistics and a figure covering the four gait parameters.

function summary = Summarise_Measured_Results(resultsFolder)
T = readtable(fullfile(resultsFolder, 'Results.csv'), 'TextType', 'string');
reference = {'gtSpeed_mps','gtStepCount','gtStrideDuration_s','gtStrideLength_m'};
estimate = {'estSpeed_mps','estStepCount','estStrideDuration_mean_s','estStrideLength_mean_m'};
labels = ["Walking speed"; "Step count"; "Stride duration"; "Stride length"];
units = ["m/s"; "steps"; "s"; "m"];
N = zeros(4, 1); MAE = NaN(4, 1); RMSE = NaN(4, 1); MAPE_pct = NaN(4, 1);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 750]);
figCleanup = onCleanup(@() close(fig));
for k = 1:4
    truth = T.(reference{k});
    estimated = T.(estimate{k});
    valid = isfinite(truth) & isfinite(estimated);
    error = estimated(valid)-truth(valid);
    N(k) = sum(valid);
    MAE(k) = mean(abs(error));
    RMSE(k) = sqrt(mean(error.^2));
    percentageValid = valid & truth ~= 0;
    MAPE_pct(k) = 100*mean(abs((estimated(percentageValid)-truth(percentageValid))./truth(percentageValid)));
    subplot(2, 2, k);
    plot(1:height(T), truth, 'ko-', 1:height(T), estimated, 'rx--', 'LineWidth', 1.2);
    xlabel('Recording'); ylabel(labels(k)+' ('+units(k)+')');
    title(sprintf('RMSE %.4g %s', RMSE(k), units(k)));
    legend('Reference', 'Estimate', 'Location', 'best');
    grid on;
end
summary = table(labels, units, N, MAE, RMSE, MAPE_pct, ...
    'VariableNames', {'Parameter','Unit','N','MAE','RMSE','MAPE_pct'});
writetable(summary, fullfile(resultsFolder, 'Summary.csv'));
exportgraphics(fig, fullfile(resultsFolder, 'Parameter_Comparison.png'), 'Resolution', 150);
disp(summary);
end
