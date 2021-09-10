function plot_evolution( Grid, U, xrange )

X = Grid{1}(1, :);
T = Grid{2}(:, 1);

sel = ((X > xrange(1)) & (X < xrange(2)));

image([T(1) T(end)], [xrange(1) xrange(2)], ...
      transpose(abs(U(:, sel))), 'CDataMapping', 'scaled')
colormap('hot')

end

