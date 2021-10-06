function plot_with_stability_smart(omegas, norms, stability)

to_plot = {omegas(1), norms(1)};
prev = stability(1);

for i = 2:length(stability)
    if prev ~= stability(i)
        if prev == 1
            plot(to_plot{1}, to_plot{2}, 'LineWidth', 3, 'Color', 'black')
        else
            plot(to_plot{1}, to_plot{2}, 'LineWidth', 0.5, 'Color', 'black')
        end
        prev = stability(i);
        to_plot = {omegas(i-1), norms(i-1)};
    end
    
    to_plot{1}(end+1) = omegas(i);
    to_plot{2}(end+1) = norms(i);
end

if prev == 1
    plot(to_plot{1}, to_plot{2}, 'LineWidth', 3, 'Color', 'black')
else
    plot(to_plot{1}, to_plot{2}, 'LineWidth', 0.5, 'Color', 'black')
end

end

