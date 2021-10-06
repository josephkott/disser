function [omegas, norms, c_parameters, stability] = ...
    get_branch_with_linear_analog( ...
        params, n, c, xspan, omega_end, norm_threshold, with_stability)

if mod(n, 2) == 0
    type = 'even';
    get_end = @(params, c) get_ux_end_cosine_nho(params, xspan, c, 8192);
else
    type = 'odd';
    get_end = @(params, c) get_u_end_cosine_nho(params, xspan, c, 8192);
end

if params(3) >= 0
    omega_step = -0.005;
else
    omega_step = 0.01;
end

omega = (2 * n) + 1 + 0.01 * sign(omega_step);
params(1) = omega;
omegas = omega:omega_step:omega_end;

% fix "c"
f = @(c) get_end(params, c);
c_init = newton(f, c);

[norms, c_parameters] = numerical_continuation( ...
    type, params, xspan, omegas, c_init, norm_threshold);

stability = zeros(1, length(omegas));

if with_stability
    for i = 1:length(omegas)
        fprintf('%g/%g\n', i, length(omegas))
        params(1) = omegas(i);
        [u0, du0] = asympt_left_cosine_nho(params, c_parameters(i), xspan(1));
        [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
        X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];
        [stability(i), max_real, min_abs] = is_stable_nho(params, X, U, 128);
        fprintf('max: %g, min: %g\n', max_real, min_abs)
    end
end

end
