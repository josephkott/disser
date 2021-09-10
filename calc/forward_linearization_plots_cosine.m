function forward_linearization_plots_cosine(params, u_span, du_span, M, N)

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
Um = f_scan_cosine(params, [0, -pi], u_span, du_span, M, N, 1024);

Island = ones(N, M);
Island(Up == 0 & Um == 0) = 0;

Upb = f_scan_bounded_cosine(params, [0, pi], u_span, du_span, M, N, 4, 20, 1024);
U2p = f_scan_cosine(params, [0, 2*pi], u_span, du_span, M, N, 1024);

U2pb = ones(N, M);
U2pb(Upb == 0 & U2p == 0) = 0;

Strips = ones(N, M);
Strips(Island == 0 & U2pb == 0) = 0;

selector = find(Strips == 0);
x = linspace( u_span(1),  u_span(2), M);
y = linspace(du_span(1), du_span(2), N);
[X, Y] = meshgrid(x, y);
Xs = X(selector);
Ys = Y(selector);

a11 = zeros(1, length(selector));
Types = Strips;
Types(Types == 1) = -1;
for i = 1:length(selector)
    index = selector(i);
    
    p = [Xs(i) Ys(i)];
    J = jacobian_f_map_cosine(params, p, 1024);
    
    Types(index) = get_signs_configuration_type(J);
    a11(i) = J(1, 1);
end

plot_scan(u_span, du_span, Types, true, 'gray');
plot_hist(a11)

end

