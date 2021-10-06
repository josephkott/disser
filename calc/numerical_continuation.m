% function [norms, c_parameters, stability] = ...
function [norms, c_parameters] = ...
    numerical_continuation(type, params, xspan, omegas, c_init, norm_threshold)

fprintf('|   |-> %g/%g\n', 1, length(omegas))
c_parameters = zeros(1, length(omegas)) * NaN;
norms = zeros(1, length(omegas)) * NaN;
% stability =  zeros(1, length(omegas)) * NaN;
c_parameters(1) = c_init;

[u0, du0] = asympt_left_cosine_nho(params, c_init, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);

if strcmp(type, 'odd')
    X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];
elseif strcmp(type, 'even')
    X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];
end

if strcmp(type, 'odd')
	get_end = @(params, c) get_u_end_cosine_nho(params, xspan, c, 8192);
elseif strcmp(type, 'even')
    get_end = @(params, c) get_ux_end_cosine_nho(params, xspan, c, 8192);
end

norms(1) = trapz(X, U .^ 2);
% [stability(1), max_real, min_abs] = is_stable_nho(params, X, U, 256);
% fprintf('|   |   | max: %g, min: %g\n', max_real, min_abs)

for i = 2:length(omegas)
    fprintf('|   |-> %g/%g\n', i, length(omegas))
    params(1) = omegas(i);
    
    f = @(c) get_end(params, c);
    c_parameters(i) = newton(f, c_parameters(i - 1));
    
    if isnan(c_parameters(i))
        disp 'Oups'
        break
    end
    
    [u0, du0] = asympt_left_cosine_nho(params, c_parameters(i), xspan(1));
    [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
    
    if strcmp(type, 'odd')
        X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];
    elseif strcmp(type, 'even')
        X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];
    end
    
	norms(i) = trapz(X, U .^ 2);
%     [stability(i), max_real, min_abs] = is_stable_nho(params, X, U, 256);
%     fprintf('|   |   | max: %g, min: %g\n', max_real, min_abs)
    
	if norms(i) > norm_threshold
        disp 'Norm'
        norms(i) = NaN;
        c_parameters(i) = NaN;
        break
	else
end

end

