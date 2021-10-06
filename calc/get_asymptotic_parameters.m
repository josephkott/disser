function [c_odds, c_odd_norms, c_evens, c_even_norms] = ...
    get_asymptotic_parameters(params, xspan, c, u_values, du_values)

eps = 1e-9;
end_threshold = 50;
norm_threshold = 50;

 u_values(abs( u_values) > end_threshold) = 0;
du_values(abs(du_values) > end_threshold) = 0;

dichotomy_borders_u  = get_dichotomy_pairs(c,  u_values); % for odd solutions
dichotomy_borders_du = get_dichotomy_pairs(c, du_values); % for even solutions

% Odd solutions
c_odds = zeros(1, length(dichotomy_borders_u));

get_u_end_cosine_nho_params = @(c) get_u_end_cosine_nho(params, xspan, c, 8192);
for i = 1:length(dichotomy_borders_u)
    a = dichotomy_borders_u{i}(1);
    b = dichotomy_borders_u{i}(2);
    c_odds(i) = dichotomy(get_u_end_cosine_nho_params, a, b, eps);
end

% Even solutions
c_evens = zeros(1, length(dichotomy_borders_du));

get_ux_end_cosine_nho_params = @(c) get_ux_end_cosine_nho(params, xspan, c, 8192);
for i = 1:length(dichotomy_borders_du)
    a = dichotomy_borders_du{i}(1);
    b = dichotomy_borders_du{i}(2);
    c_evens(i) = dichotomy(get_ux_end_cosine_nho_params, a, b, eps);
end

% Norms for odd solutions
c_odd_norms = zeros(1, length(c_odds));
for i = 1:length(c_odds)
    [u0, du0] = asympt_left_cosine_nho(params, c_odds(i), xspan(1));
    [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
    X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];
    c_odd_norms(i) = trapz(X, U .^ 2);
end

% Norms for even solutions
c_even_norms = zeros(1, length(c_evens));
for i = 1:length(c_evens)
    [u0, du0] = asympt_left_cosine_nho(params, c_evens(i), xspan(1));
    [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
    X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];
    c_even_norms(i) = trapz(X, U .^ 2);
end

sel = c_odd_norms  < norm_threshold;
c_odds = c_odds(sel);
c_odd_norms = c_odd_norms(sel);

sel = c_even_norms  < norm_threshold;
c_evens = c_evens(sel);
c_even_norms = c_even_norms(sel);

end

