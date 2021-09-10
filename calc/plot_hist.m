function plot_hist(values)

figure('Position', [100, 100, 200, 120])
histogram(log(abs(values)), 'Normalization', 'probability', 'NumBins', 20)
min_value_log = min(log(abs(values)));
min_value     = min(abs(values));

xticks([0 min_value_log])
yticks([])

xticklabels({'0', min_value})
xlim([0 max(log(abs(values)))])

end

