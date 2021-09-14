function plot_with_stability(omegas, norms, stability)

for i = 2:length(stability)
    if stability(i) == 1
        plot([omegas(i-1) omegas(i)], [norms(i-1), norms(i)], ...
             'LineWidth', 3, 'Color', 'black')
    else
        plot([omegas(i-1) omegas(i)], [norms(i-1), norms(i)], ...
             'LineWidth', 0.5, 'Color', 'black')
    end
end

end

