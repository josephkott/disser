function backward_linearization_plots_piecewise(params, u_span, du_span, M, N)

L1 = params(2); L2 = params(3);
L = L1 + L2;

Up = f_scan_piecewise(params, [0, L1], u_span, du_span, M, N, 1024);
Um = f_scan_piecewise(params, [0, -L], u_span, du_span, M, N, 1024);

Island = ones(N, M);
Island(Up == 0 & Um == 0) = 0;

Umb = f_scan_bounded_piecewise(params, [0, -L], u_span, du_span, M, N, 2, 4, 1024);
U2m = f_scan_piecewise(params, [0, -2*L], u_span, du_span, M, N, 1024);

U2mb = ones(N, M);
U2mb(Umb == 0 & U2m == 0) = 0;

Strips = ones(N, M);
Strips(Island == 0 & U2mb == 0) = 0;

selector = find(Strips == 0);
x = linspace( u_span(1),  u_span(2), M);
y = linspace(du_span(1), du_span(2), N);
[X, Y] = meshgrid(x, y);
Xs = X(selector);
Ys = Y(selector);

b22 = zeros(1, length(selector));
Types = Strips;
Types(Types == 1) = -1;
for i = 1:length(selector)
    index = selector(i);
    
    p = [Xs(i) Ys(i)];
    J = jacobian_f_imap_piecewise(params, p, 1024);
    
    Types(index) = get_signs_configuration_type(J);
    b22(i) = J(2, 2);
end

plot_scan(u_span, du_span, Types, true, 'gray');
plot_hist(b22)

end
