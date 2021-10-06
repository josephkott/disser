%% Picture (1). Localized solutions with and without linear counterpart, $P_0 \neq 0$
clc; clear

addpath('./mex/')

omega = 8; T = 12; sigma0 = 1; P1 = 2;
params = [omega T sigma0 P1];
xspan = [-6 0];

c = 0:0.5:250;
left = zeros(length(c), 2);

for i = 1:length(c)
	[u0, du0] = asympt_left_cosine_nho(params, c(i), xspan(1));
	[X, U] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
	left(i, :) = U(end, :);
end

figure; plot(c, left(:, 1))

%% Collect localized solutions for Picture (1).

eps = 1e-9;
end_threshold = 50;
norm_threshold = 50;

u_values  = left(:, 1);  u_values(abs( u_values) > end_threshold) = 0;
du_values = left(:, 2); du_values(abs(du_values) > end_threshold) = 0;

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

c_odds  = c_odds (c_odd_norms  < norm_threshold);
c_evens = c_evens(c_even_norms < norm_threshold);

%% Numerical continuation

omegas = omega:0.01:8;

% Odd solutions
fprintf('odd solutions:\n')
norms_odd = zeros(length(c_odds), length(omegas));
c_parameters_odd = zeros(length(c_odds), length(omegas));

for i = 1:length(c_odds)
    fprintf('\t%g/%g\n', i, length(c_odds))
    [norms, c_parameters] = ...
        numerical_continuation('odd', params, xspan, omegas, c_odds(i), norm_threshold);
    
    norms_odd(i, :) = norms;
    c_parameters_odd(i, :) = c_parameters;
end

% Even solutions
fprintf('even solutions:\n')
norms_even = zeros(length(c_evens), length(omegas));
c_parameters_even = zeros(length(c_evens), length(omegas));

for i = 1:length(c_evens)
    fprintf('\t%g/%g\n', i, length(c_evens))
    [norms, c_parameters] = ...
        numerical_continuation('even', params, xspan, omegas, c_evens(i), norm_threshold);
    
    norms_even(i, :) = norms;
    c_parameters_even(i, :) = c_parameters;
end

%% Plot them all!

figure; hold on

for i = 1:length(c_odds)
    plot(omegas, norms_odd(i, :))
end

for i = 1:length(c_evens)
    plot(omegas, norms_even(i, :))
end

%% Moving from the left to right

clc; clear

addpath('./mex/')

T = 12; sigma0 = 1; P1 = 2;
omega_1st_steps = -2:0.5:8;

xspan = [-6 0];
c = 0:0.5:250;

eps = 1e-9;

figure; hold on % tracker
parts = {};
for i = 1:(length(omega_1st_steps) - 1)
    fprintf('omega = %g (%g/%g)\n', omega_1st_steps(i), i, length(omega_1st_steps))
    omega = omega_1st_steps(i);
    params = [omega T sigma0 P1];
    
    left = zeros(length(c), 2);
    for j = 1:length(c)
        [u0, du0] = asympt_left_cosine_nho(params, c(j), xspan(1));
        [~, U] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
        left(j, :) = U(end, :);
    end
    
    % TODO: add norm_threshold to parameters
    [c_odds, c_odd_norms, c_evens, c_even_norms] = ...
        get_asymptotic_parameters(params, xspan, c, left(:, 1), left(:, 2));
    
    omegas = omega_1st_steps(i):0.0025:omega_1st_steps(i + 1);
    
    % Odd solutions continuation
    fprintf('|-> odd\n')
    for j = 1:length(c_odds)
        fprintf('|   |-> %g/%g\n', j, length(c_odds))
        [norms, c_parameters] = numerical_continuation( ...
            'odd', params, xspan, omegas, c_odds(j), norm_threshold);
        parts{end + 1} = {omegas, norms, c_parameters};
        pause(1e-6); plot(omegas, norms, 'Color', 'black'); pause(1e-6) % tracker
    end
    
    % Even solutions continuation
    fprintf('|-> even\n')
    for j = 1:length(c_evens)
        fprintf('|   |-> %g/%g\n', j, length(c_evens))
        [norms, c_parameters] = numerical_continuation( ...
            'even', params, xspan, omegas, c_evens(j), norm_threshold);
        parts{end + 1} = {omegas, norms, c_parameters};
        pause(1e-6); plot(omegas, norms, 'Color', 'black'); pause(1e-6) % tracker
    end
end

%%

figure; hold on
for i = 1:length(parts)
    plot3(parts{i}{1}, parts{i}{2}, ones(1, length(parts{i}{1})) * i, 'Color', 'black');
end
view(2)

%% Picture (1). Without linear counterpart
omega_1a = 5; c_parameter_1a = 53.693031313829124; type_1a = 'odd';
omega_1b = 5; c_parameter_1b = 46.753605435602367; type_1b = 'odd';

omega_2a = 3.5; c_parameter_2a = 1.205846722284332e+02; type_2a = 'even';
omega_2b = 3; c_parameter_2b = 1.785822237106040e+02; type_2b = 'even';

omega_3a = -0.5; c_parameter_3a = 49.431950644589961; type_3a = 'odd';
omega_3b = -0.5; c_parameter_3b = 21.234334138222039; type_3b = 'odd';

omega_4a = -1; c_parameter_4a = 68.606987769715488; type_4a = 'even';
omega_4b = -1.5; c_parameter_4b = 96.971709842793643; type_4b = 'even';

omegas_start = [
    omega_1a omega_1b, ...
    omega_2a omega_2b, ...
    omega_3a omega_3b, ...
    omega_4a omega_4b, ...
];

c_parameters_start = [
    c_parameter_1a c_parameter_1b, ...
    c_parameter_2a c_parameter_2b, ...
    c_parameter_3a c_parameter_3b, ...
    c_parameter_4a c_parameter_4b, ...
];

types = {
    type_1a type_1b, ...
    type_2a type_2b, ...
    type_3a type_3b, ...
    type_4a type_4b, ...
};

collection = {};
omega_step = 0.0025;
norm_threshold = 30;
for i = 1:8
    fprintf('%g/%g\n', i, 8)
    
    % to the right
    fprintf('|-> right\n')
    omegas = omegas_start(i):omega_step:8;
    params(1) = omegas_start(i);
    [norms, c_parameters] = numerical_continuation( ...
        types{i}, params, xspan, omegas, c_parameters_start(i), norm_threshold);
    
    collection{end + 1} = {omegas, norms, c_parameters};

    % to the left
    fprintf('|-> left\n')
    omegas = omegas_start(i):(-omega_step):-2;
    params(1) = omegas_start(i);
    [norms, c_parameters] = numerical_continuation( ...
        types{i}, params, xspan, omegas, c_parameters_start(i), norm_threshold);
    
    collection{end + 1} = {omegas, norms, c_parameters};
end

figure('Position', [100, 100, 400, 300]); hold on
for i = 1:length(collection)
    plot(collection{i}{1}, collection{i}{2}, 'Color', 'black')
end

axis([-2 8 0 22])
xticks(-2:1:8)

%% Solutions without linear counterpart (a)

params(1) = collection{15}{1}(1); % -1.5
c_parameter = collection{15}{3}(1); % 96.971709842793643
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (a)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'black')

%% Solutions without linear counterpart (b)

params(1) = collection{9}{1}(1001); % 0.5
c_parameter = collection{9}{3}(1001); % 29.977276881176259
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (b)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'black')

%% Solutions without linear counterpart (c)

params(1) = collection{8}{1}(901); % 2.1
c_parameter = collection{8}{3}(901); % 3.122441206154584e+02
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (c)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'black')

%% Solutions without linear counterpart (d)

params(1) = collection{2}{1}(501); % 4.5
c_parameter = collection{2}{3}(501); % 72.880679079470951
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (d)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'black')

%% Picture (1). With linear counterpart

figure('Position', [100, 100, 400, 300]); hold on
axis([-2 8 0 22])
xticks(-2:1:8)

[omegas_G0, norms_G0, c_parameters_G0, stability_G0] = ...
    get_branch_with_linear_analog(params, 0, 0.5, xspan, -2, false);
% plot_with_stability(omegas_G0, norms_G0, stability_G0)
plot(omegas_G0, norms_G0, 'Color', 'black')

[omegas_G1, norms_G1, c_parameters_G1, stability_G1] = ...
    get_branch_with_linear_analog(params, 1, 0.5, xspan, -2, false);
% plot_with_stability(omegas_G1, norms_G1, stability_G1)
plot(omegas_G1, norms_G1, 'Color', 'black')

[omegas_G2, norms_G2, c_parameters_G2, stability_G2] = ...
    get_branch_with_linear_analog(params, 2, 0.5, xspan, -2, false);
% plot_with_stability(omegas_G2, norms_G2, stability_G2)
plot(omegas_G2, norms_G2, 'Color', 'black')

[omegas_G3, norms_G3, c_parameters_G3, stability_G3] = ...
    get_branch_with_linear_analog(params, 3, 0.5, xspan, -2, false);
plot(omegas_G3, norms_G3, 'Color', 'black')
% plot_with_stability(omegas_G3, norms_G3, stability_G3)

%% Solution with linear counterpart, $\Gamma_0$

params(1) = omegas_G0(100);
c_parameter = c_parameters_G0(100);
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -0.2 1.4])

%% Set tip, $\Gamma_0$
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solution with linear counterpart, $\Gamma_1$

params(1) = omegas_G1(100);
c_parameter = c_parameters_G1(100);
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip, $\Gamma_1$
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solution with linear counterpart, $\Gamma_2$

params(1) = omegas_G2(100);
c_parameter = c_parameters_G2(100);
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip, $\Gamma_2$
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solution with linear counterpart, $\Gamma_3$

params(1) = omegas_G3(100);
c_parameter = c_parameters_G3(100);
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip, $\Gamma_3$
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Picture (1b). P0 < 0

omega = 8; T = 12; sigma0 = -1; P1 = 2;
params = [omega T sigma0 P1];
xspan = [-6 0];

c = 0:0.02:25;
left = zeros(length(c), 2);

for i = 1:length(c)
	[u0, du0] = asympt_left_cosine_nho(params, c(i), xspan(1));
	[X, U] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
	left(i, :) = U(end, :);
end

figure; plot(c, left(:, 1))

%% Picture (1b).

clc; clear

addpath('./mex/')

omega = 8; T = 12; sigma0 = -1; P1 = 2;
params = [omega T sigma0 P1];
xspan = [-6 0];

c = 0:0.01:20;
left = zeros(length(c), 2);

for i = 1:length(c)
	[u0, du0] = asympt_left_cosine_nho(params, c(i), xspan(1));
	[X, U] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
	left(i, :) = U(end, :);
end

figure; plot(c, left(:, 1))

%% Moving from the right to left

clc; clear

addpath('./mex/')

T = 12; sigma0 = -1; P1 = 2;
omega_1st_steps = 8:-0.5:-2;

xspan = [-6 0];
c = 0:0.01:20;

eps = 1e-9;
norm_threshold = 50;

figure; hold on % tracker
parts = {};
for i = 1:(length(omega_1st_steps) - 1)
    fprintf('omega = %g (%g/%g)\n', omega_1st_steps(i), i, length(omega_1st_steps))
    omega = omega_1st_steps(i);
    params = [omega T sigma0 P1];
    
    left = zeros(length(c), 2);
    for j = 1:length(c)
        [u0, du0] = asympt_left_cosine_nho(params, c(j), xspan(1));
        [~, U] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
        left(j, :) = U(end, :);
    end
    
    % TODO: add norm_threshold to parameters
    [c_odds, c_odd_norms, c_evens, c_even_norms] = ...
        get_asymptotic_parameters(params, xspan, c, left(:, 1), left(:, 2));
    
    omegas = omega_1st_steps(i):(-0.0025):omega_1st_steps(i + 1);
    
    % Odd solutions continuation
    fprintf('|-> odd\n')
    for j = 1:length(c_odds)
        fprintf('|   |-> %g/%g\n', j, length(c_odds))
        [norms, c_parameters] = numerical_continuation( ...
            'odd', params, xspan, omegas, c_odds(j), norm_threshold);
        parts{end + 1} = {omegas, norms, c_parameters};
        pause(1e-6); plot(omegas, norms, 'Color', 'black'); pause(1e-6) % tracker
    end
    
    % Even solutions continuation
    fprintf('|-> even\n')
    for j = 1:length(c_evens)
        fprintf('|   |-> %g/%g\n', j, length(c_evens))
        [norms, c_parameters] = numerical_continuation( ...
            'even', params, xspan, omegas, c_evens(j), norm_threshold);
        parts{end + 1} = {omegas, norms, c_parameters};
        pause(1e-6); plot(omegas, norms, 'Color', 'black'); pause(1e-6) % tracker
    end
end

%% Picture (1b). Without linear counterpart
omega_1 = 3;   c_parameter_1 = 2.044801527857780; type_1 = 'even';
omega_2 = 1.5; c_parameter_2 = 2.827211639285087; type_2 = 'odd';

omegas_start = [omega_1 omega_2];
c_parameters_start = [c_parameter_1 c_parameter_2];
types = {type_1 type_2};

collection = {};
omega_step = 0.01;
norm_threshold = 50;
for i = 1:2
    fprintf('%g/%g\n', i, 2)
    
    % to the right
    fprintf('|-> right\n')
    omegas = omegas_start(i):omega_step:8;
    params(1) = omegas_start(i);
    [norms, c_parameters, stability] = numerical_continuation( ...
        types{i}, params, xspan, omegas, c_parameters_start(i), norm_threshold);
    
    collection{end + 1} = {omegas, norms, c_parameters, stability};

    % to the left
    fprintf('|-> left\n')
    omegas = omegas_start(i):(-omega_step):-2;
    params(1) = omegas_start(i);
    [norms, c_parameters] = numerical_continuation( ...
        types{i}, params, xspan, omegas, c_parameters_start(i), norm_threshold);
    
    collection{end + 1} = {omegas, norms, c_parameters, stability};
end

figure('Position', [100, 100, 400, 300]); hold on
for i = 1:length(collection)
    % plot_with_stability(collection{i}{1}, collection{i}{2}, collection{i}{4})
    plot(collection{i}{1}, collection{i}{2}, 'Color', 'black')
end

axis([-2 8 0 50])
xticks(-2:1:8)

%% Solutions without linear counterpart (a)

params(1) = 3;
c_parameter = 2.044801527857780;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -1 9])

%% Set tip (a)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solutions without linear counterpart (b)

params(1) = 1.5;
c_parameter = 2.827211639285087;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (b)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Picture (1b). With linear counterpart

figure('Position', [100, 100, 400, 300]); hold on
axis([-2 8 0 50])
xticks(-2:1:8)

[omegas_G0, norms_G0, c_parameters_G0, stability_G0] = ...
    get_branch_with_linear_analog(params, 0, 0.5, xspan, 8, norm_threshold, false);
% plot_with_stability(omegas_G0, norms_G0, stability_G0)
plot(omegas_G0, norms_G0, 'Color', 'black')

[omegas_G1, norms_G1, c_parameters_G1, stability_G1] = ...
    get_branch_with_linear_analog(params, 1, 0.5, xspan, 8, norm_threshold, false);
% plot_with_stability(omegas_G1, norms_G1, stability_G1)
plot(omegas_G1, norms_G1, 'Color', 'black')

[omegas_G2, norms_G2, c_parameters_G2, stability_G2] = ...
    get_branch_with_linear_analog(params, 2, 0.5, xspan, 8, norm_threshold, false);
% plot_with_stability(omegas_G2, norms_G2, stability_G2)
plot(omegas_G2, norms_G2, 'Color', 'black')

[omegas_G3, norms_G3, c_parameters_G3, stability_G3] = ...
    get_branch_with_linear_analog(params, 3, 0.5, xspan, 8, norm_threshold, false);
plot(omegas_G3, norms_G3, 'Color', 'black')
% plot_with_stability(omegas_G3, norms_G3, stability_G3)

%% Picture (1b). Stability for solutions with linear cunterpart

sel = 1:10:length(c_parameters_G3);
sub_omegas = omegas_G3(sel);
sub_c_parameters = c_parameters_G3(sel);
sub_norms = norms_G3(sel);

sel = ~isnan(sub_c_parameters);
sub_omegas = sub_omegas(sel);
sub_c_parameters = sub_c_parameters(sel);
sub_norms = sub_norms(sel);

stability_G3 = zeros(1, length(sub_c_parameters));
params(4) = 2;
for i = 1:length(sub_c_parameters)
    fprintf('%g/%g\n', i, length(sub_c_parameters))
    params(1) = sub_omegas(i);
    [u0, du0] = asympt_left_cosine_nho(params, sub_c_parameters(i), xspan(1));
    [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
    X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];
    [stability_G3(i), max_real, min_abs] = is_stable_nho(params, X, U, 128);
    fprintf('max: %g, min: %g\n', max_real, min_abs)
end

sub_norms_G3 = sub_norms;
sub_omegas_G3 = sub_omegas;

figure; hold on
plot_with_stability(sub_omegas_G3, sub_norms_G3, stability_G3)

figure('Position', [100, 100, 400, 300]); hold on
plot_with_stability(sub_omegas_G0, sub_norms_G0, stability_G0)
plot_with_stability(sub_omegas_G1, sub_norms_G1, stability_G1)
plot_with_stability(sub_omegas_G2, sub_norms_G2, stability_G2)
plot_with_stability(sub_omegas_G3, sub_norms_G3, stability_G3)
axis([-2 8 0 50])
xticks(-2:1:8)

%% Solution $\Gamma_0$

index = find(omegas_G0 == 3);
params(1) = omegas_G0(index);
c_parameter = c_parameters_G0(index);

[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -0.5 2])

%% Solution $\Gamma_1$

index = find(omegas_G1 == 5);
params(1) = omegas_G1(index);
c_parameter = c_parameters_G1(index);

[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Solution $\Gamma_2$

index = find(omegas_G2 == 7);
params(1) = omegas_G2(index);
c_parameter = c_parameters_G2(index);

[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Solution $\Gamma_3$

index = find(omegas_G3 == 8);
params(1) = omegas_G3(index);
c_parameter = c_parameters_G3(index);

[u0, du0] = asympt_left_cosine_nho(params, c_parameter, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Picture (2)

clc; clear

addpath('./mex/')

omega = -2; T = 8; sigma0 = 0; P1 = 1;
params = [omega T sigma0 P1];
xspan = [-6 0];

c = 0:0.01:20;
left = zeros(length(c), 2);

for i = 1:length(c)
	[u0, du0] = asympt_left_cosine_nho(params, c(i), xspan(1));
	[X, U] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
	left(i, :) = U(end, :);
end

figure; plot(c, left(:, 1))

%% Moving from the left to right

clc; clear

addpath('./mex/')

T = 8; sigma0 = 0; P1 = 1;
omega_1st_steps = -4:0.5:0;

xspan = [-6 0];
c = 0:0.01:100;

eps = 1e-9;
norm_threshold = 50;

figure; hold on % tracker
parts = {};
for i = 1:(length(omega_1st_steps) - 1)
    fprintf('omega = %g (%g/%g)\n', omega_1st_steps(i), i, length(omega_1st_steps))
    omega = omega_1st_steps(i);
    params = [omega T sigma0 P1];
    
    left = zeros(length(c), 2);
    for j = 1:length(c)
        [u0, du0] = asympt_left_cosine_nho(params, c(j), xspan(1));
        [~, U] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
        left(j, :) = U(end, :);
    end
    
    % TODO: add norm_threshold to parameters
    [c_odds, c_odd_norms, c_evens, c_even_norms] = ...
        get_asymptotic_parameters(params, xspan, c, left(:, 1), left(:, 2));
    
    omegas = omega_1st_steps(i):0.0025:omega_1st_steps(i + 1);
    
    % Odd solutions continuation
    fprintf('|-> odd\n')
    for j = 1:length(c_odds)
        fprintf('|   |-> %g/%g\n', j, length(c_odds))
        [norms, c_parameters] = numerical_continuation( ...
            'odd', params, xspan, omegas, c_odds(j), norm_threshold);
        parts{end + 1} = {omegas, norms, c_parameters};
        pause(1e-6); plot(omegas, norms, 'Color', 'black'); pause(1e-6) % tracker
    end
    
    % Even solutions continuation
    fprintf('|-> even\n')
    for j = 1:length(c_evens)
        fprintf('|   |-> %g/%g\n', j, length(c_evens))
        [norms, c_parameters] = numerical_continuation( ...
            'even', params, xspan, omegas, c_evens(j), norm_threshold);
        parts{end + 1} = {omegas, norms, c_parameters};
        pause(1e-6); plot(omegas, norms, 'Color', 'black'); pause(1e-6) % tracker
    end
end

%%

figure; hold on
for i = 1:length(parts)
    plot3(parts{i}{1}, parts{i}{2}, ones(1, length(parts{i}{1})) * i, 'Color', 'black');
end
view(2)

%% Picture (1). Without linear counterpart
omega_1a = 3.5; c_parameter_1a = 8.021231346726417; type_1a = 'odd';
omega_1b = 4; c_parameter_1b = 16.646330961585051; type_1b = 'odd';

omega_2a = 0.5; c_parameter_2a = 14.607061173319813; type_2a = 'even';
omega_2b = 0.5; c_parameter_2b = 6.208928976655006; type_2b = 'even';

omega_3a = -0.5; c_parameter_3a = 24.133432341217997; type_3a = 'odd';
omega_3b = -0.5; c_parameter_3b = 22.905353004336362; type_3b = 'odd';

omega_4a = 5; c_parameter_4a = 14.230072169899941; type_4a = 'even';
omega_4b = 5.5; c_parameter_4b = 11.286635268330578; type_4b = 'even';

omegas_start = [
    omega_1a omega_1b, ...
    omega_2a omega_2b, ...
    omega_3a omega_3b, ...
    omega_4a omega_4b, ...
];

c_parameters_start = [
    c_parameter_1a c_parameter_1b, ...
    c_parameter_2a c_parameter_2b, ...
    c_parameter_3a c_parameter_3b, ...
    c_parameter_4a c_parameter_4b, ...
];

types = {
    type_1a type_1b, ...
    type_2a type_2b, ...
    type_3a type_3b, ...
    type_4a type_4b, ...
};

collection = {};
omega_step = 0.0025;
% omega_step = 0.001;
norm_threshold = 50;
for i = 1:length(omegas_start)
    fprintf('%g/%g\n', i, length(omegas_start))
    
    % to the right
    fprintf('|-> right\n')
    omegas = omegas_start(i):omega_step:8;
    params(1) = omegas_start(i);
    [norms, c_parameters] = numerical_continuation( ...
        types{i}, params, xspan, omegas, c_parameters_start(i), norm_threshold);
    
    collection{end + 1} = {omegas, norms, c_parameters};

    % to the left
    fprintf('|-> left\n')
    omegas = omegas_start(i):(-omega_step):-2;
    params(1) = omegas_start(i);
    [norms, c_parameters] = numerical_continuation( ...
        types{i}, params, xspan, omegas, c_parameters_start(i), norm_threshold);
    
    collection{end + 1} = {omegas, norms, c_parameters};
end

figure('Position', [100, 100, 400, 300]); hold on
for i = 1:length(collection)
    plot(collection{i}{1}, collection{i}{2}, 'Color', 'black')
end

axis([-2 8 0 50])
xticks(-2:1:8)

%% Check stability

reduced = {};
% for i = 1:length(collection)
for i = 2:2
    sel = 1:10:length(collection{i}{1});
    reduced{i} = {collection{i}{1}(sel), collection{i}{2}(sel), collection{i}{3}(sel)};

    sel = ~isnan(reduced{i}{3});
    reduced{i}{1} = reduced{i}{1}(sel);
    reduced{i}{2} = reduced{i}{2}(sel);
    reduced{i}{3} = reduced{i}{3}(sel);
    
    stability = zeros(1, length(reduced{i}{1}));
    for j = 1:length(reduced{i}{3})
        fprintf('%g/%g\n', j, length(reduced{i}{3}))
        params(1) = reduced{i}{1}(j);
        [u0, du0] = asympt_left_cosine_nho(params, reduced{i}{3}(j), xspan(1));
        [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
        if strcmp(types(floor(i / 2) + mod(i, 2)), 'even')
            X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];
        else
            X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];
        end
        [stability(j), max_real, min_abs] = is_stable_nho(params, X, U, 128);
        fprintf('max: %g, min: %g\n', max_real, min_abs)
    end
    
    reduced{i}{4} = stability;
end

%%

figure('Position', [100, 100, 400, 300]); hold on
for i = 1:16
    plot_with_stability(reduced{i}{1}, reduced{i}{2}, reduced{i}{4})
end

%% Solution (a)

params(1) = omega_3b;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter_3b, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -5.5 5.5])

%% Set tip (a)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solution (a')

params(1) = omega_3a;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter_3a, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -5.5 5.5])

%% Set tip (a')
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solution (b)

params(1) = omega_2a;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter_2a, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (b)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solution (b')

params(1) = omega_2b;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter_2b, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (b')
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Solution (c)

params(1) = omega_1a;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter_1a, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on

%% Set tip (c)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Soltuion (d)

params(1) = omega_4b;
[u0, du0] = asympt_left_cosine_nho(params, c_parameter_4b, xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -5.5 5.5])

%% Set tip (d)
norm = trapz(X, U .^ 2);
plot(params(1), norm, '.', 'MarkerSize', 10, 'Color', 'red')

%% Picture (2). With linear counterpart

figure('Position', [100, 100, 400, 300]); hold on
axis([-2 8 0 50])
xticks(-2:1:8)

[omegas_G0, norms_G0, c_parameters_G0, stability_G0] = ...
    get_branch_with_linear_analog(params, 0, 1, xspan, 0, norm_threshold, false);
% plot_with_stability(omegas_G0, norms_G0, stability_G0)
plot(omegas_G0, norms_G0, 'Color', 'black')

[omegas_G1, norms_G1, c_parameters_G1, stability_G1] = ...
    get_branch_with_linear_analog(params, 1, 1, xspan, 2, norm_threshold, false);
% plot_with_stability(omegas_G1, norms_G1, stability_G1)
plot(omegas_G1, norms_G1, 'Color', 'black')

[omegas_G2, norms_G2, c_parameters_G2, stability_G2] = ...
    get_branch_with_linear_analog(params, 2, 1, xspan, 4, norm_threshold, false);
% plot_with_stability(omegas_G2, norms_G2, stability_G2)
plot(omegas_G2, norms_G2, 'Color', 'black')

[omegas_G3, norms_G3, c_parameters_G3, stability_G3] = ...
    get_branch_with_linear_analog(params, 3, 2, xspan, 6, norm_threshold, false);
% plot_with_stability(omegas_G3, norms_G3, stability_G3)
plot(omegas_G3, norms_G3, 'Color', 'black')

%%

sel = 1:50:length(c_parameters_G3);
sub_omegas = omegas_G3(sel);
sub_c_parameters = c_parameters_G3(sel);
sub_norms = norms_G3(sel);

sel = ~isnan(sub_c_parameters);
sub_omegas = sub_omegas(sel);
sub_c_parameters = sub_c_parameters(sel);
sub_norms = sub_norms(sel);

sub_stability = zeros(1, length(sub_c_parameters));
for i = 1:length(sub_c_parameters)
    fprintf('%g/%g\n', i, length(sub_c_parameters))
    params(1) = sub_omegas(i);
    [u0, du0] = asympt_left_cosine_nho(params, sub_c_parameters(i), xspan(1));
    [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
    X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];
    [sub_stability(i), max_real, min_abs] = is_stable_nho(params, X, U, 128);
    fprintf('max: %g, min: %g\n', max_real, min_abs)
end

%% Solution with linear counterpart (a)

params(1) = omegas_G0(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G0(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -0.5 4])

%% Solution with linear counterpart (b)

params(1) = omegas_G1(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G1(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -3.5 3.5])

%% Solution with linear counterpart (c)

params(1) = omegas_G2(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G2(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -3 3.5])

%% Solution with linear counterpart (d)

params(1) = omegas_G3(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G3(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh(end:-1:1)]; U = [Uh(:, 1); -Uh(end:-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -3.5 3.5])

%% Picture (3)
clc; clear

addpath('./mex/')

omega = 8; T = 12; sigma0 = 1; P1 = 2;
params = [omega T sigma0 P1];
xspan = [-6 0];

%% Picutre (3), stable solution (G1, omega = 2)

[omegas_G1, norms_G1, c_parameters_G1, stability_G1] = ...
    get_branch_with_linear_analog(params, 1, 0.5, xspan, 2, 50, false);

params(1) = omegas_G1(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G1(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); -Uh((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -1.5 1.5])

eigenvalues = get_spectrum_nho(params, X, U, 256);
figure('Position', [100, 100, 350, 250])
plot_spectrum_nho(eigenvalues)

[Grid, Solution, ~] = CFDS_nho(params, X, U, 80, 0.01);
figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-6 6])

%% Picture (3), unstable solution (G1, omega = 0)

[omegas_G1, norms_G1, c_parameters_G1, stability_G1] = ...
    get_branch_with_linear_analog(params, 1, 0.5, xspan, 0, 50, false);

params(1) = omegas_G1(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G1(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); -Uh((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -2.5 2.5])

eigenvalues = get_spectrum_nho(params, X, U, 256);
figure('Position', [100, 100, 350, 250])
axis([-0.3 0.3 -3 3])
plot_spectrum_nho(eigenvalues)

[Grid, Solution, Norm] = CFDS_nho(params, X, U, 80, 0.01);
figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-6 6])

%% Picture (4)
clc; clear

addpath('./mex/')

omega = -2; T = 8; sigma0 = 0; P1 = 1;
params = [omega T sigma0 P1];

xspan = [-6 0];

%% Picutre (4), stable solution (G1, omega = 2)

[omegas_G1, norms_G1, c_parameters_G1, stability_G1] = ...
    get_branch_with_linear_analog(params, 1, 2, xspan, 2, 50, false);

params(1) = omegas_G1(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G1(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); -Uh((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -3.5 3.5])

eigenvalues = get_spectrum_nho(params, X, U, 256);
figure('Position', [100, 100, 350, 250])
plot_spectrum_nho(eigenvalues)

[Grid, Solution, ~] = CFDS_nho(params, X, U, 600, 0.05);
figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-6 6])

%% Picture (4), unstable solution (G1, omega = 2.5)

[omegas_G1, norms_G1, c_parameters_G1, stability_G1] = ...
    get_branch_with_linear_analog(params, 1, 1, xspan, 2.5, 50, false);

params(1) = omegas_G1(end);
[u0, du0] = asympt_left_cosine_nho(params, c_parameters_G1(end), xspan(1));
[Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); -Uh((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')
grid on
axis([-6 6 -3 3])

eigenvalues = get_spectrum_nho(params, X, U, 256);
figure('Position', [100, 100, 350, 250])
axis([-0.3 0.3 -3 3])
plot_spectrum_nho(eigenvalues)
axis([-0.3 0.3 -15 15])

[Grid, Solution, ~] = CFDS_nho(params, X, U, 600, 0.05);
figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-6 6])

%% Picture (5), branches with linear counterpart for different T values, non-zer mean nho
clc; clear

addpath('./mex/')

omega = 8; T = [8 12 16 0]; sigma0 = 1; P1 = [2 2 2 0];
xspan = [-6 0];

collection = {};
for i = 1:length(T)
    params = [omega T(i) sigma0 P1(i)];
    
    [omegas, norms, c_parameters, ~] = ...
        get_branch_with_linear_analog(params, 2, 1, xspan, 2.5, 50, false);
    collection{i} = {omegas, norms, c_parameters};
end

figure('Position', [100, 100, 400, 300]); hold on
axis([-2 8 0 6])
xticks(-2:1:8)

for i = 1:length(collection)
    plot(collection{i}{1}, collection{i}{2}, 'Color', 'black');
end

% Add stability
reduced = {};

type = 'even';
for i = 1:length(collection)
    fprintf('Branch: %g/%g\n', i, length(collection))
    sel = 1:10:length(collection{i}{1});
    
    omegas = collection{i}{1}(sel);
    norms = collection{i}{2}(sel);
    c_parameters = collection{i}{3}(sel);
    stability = zeros(1, length(sel));
    
    params = [omega T(i) sigma0 P1(i)];
    for j = 1:length(sel)
        fprintf('--> %g/%g\n', j, length(sel))
        params(1) = omegas(j);
        [u0, du0] = asympt_left_cosine_nho(params, c_parameters(j), xspan(1));
        [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
        if strcmp(type, 'even')
            X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); Uh((end-1):-1:1, 1)];
        else
            X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); -Uh((end-1):-1:1, 1)];
        end
        [stability(j), max_real, min_abs] = is_stable_nho(params, X, U, 128);
        fprintf('max: %g, min: %g\n', max_real, min_abs)
    end
    
    reduced{i} = {omegas, norms, c_parameters, stability};
end

reduced_G2 = reduced;

figure('Position', [100, 100, 400, 200]); hold on
axis([-2 6 0 8])
xticks(-2:1:6)

for i = 1:length(reduced_G0)
    plot_with_stability_smart(reduced_G0{i}{1}, reduced_G0{i}{2}, reduced_G0{i}{4});
end

for i = 1:length(reduced_G1)
    plot_with_stability_smart(reduced_G1{i}{1}, reduced_G1{i}{2}, reduced_G1{i}{4});
end

for i = 1:length(reduced_G2)
    plot_with_stability_smart(reduced_G2{i}{1}, reduced_G2{i}{2}, reduced_G2{i}{4});
end

figure('Position', [100, 100, 400, 200]); hold on
axis([-2 6 0 8])
xticks(-2:1:6)

for i = 1:length(reduced_G0)
    plot_with_stability_smart(reduced_G0{i}{1}, i * ones(1, length(reduced_G0{i}{2})), reduced_G0{i}{4});
end

for i = 1:length(reduced_G1)
    sel = reduced_G1{i}{1} > 1.5;
    plot_with_stability_smart(reduced_G1{i}{1}(sel), i * ones(1, length(reduced_G1{i}{2}(sel))), reduced_G1{i}{4}(sel));
end

for i = 1:length(reduced_G2)
    sel = reduced_G2{i}{1} > 3.5;
    plot_with_stability_smart(reduced_G2{i}{1}(sel), i * ones(1, length(reduced_G2{i}{2}(sel))), reduced_G2{i}{4}(sel));
end

%% Picture (6), branches with linear counterpart for different T values, non-zer mean nho

clc; clear

addpath('./mex/')

omega = 8; T = [8 12 16]; P0 = 0; sigma1 = 1;
xspan = [-6 0];

collection = {};
for i = 1:length(T)
    params = [omega T(i) P0 sigma1];
    
    [omegas, norms, c_parameters, ~] = ...
        get_branch_with_linear_analog(params, 3, 2, xspan, 5, 70, false);
    collection{i} = {omegas, norms, c_parameters};
end

figure('Position', [100, 100, 400, 200]); hold on
axis([-1 8 0 60])
xticks(-1:1:8)

for i = 1:length(collection)
    plot(collection{i}{1}, collection{i}{2}, 'Color', 'black');
end

% Add stability
reduced = {};

type = 'odd';
for i = 1:length(collection)
    fprintf('Branch: %g/%g\n', i, length(collection))
    sel = 1:2:length(collection{i}{1});
    
    omegas = collection{i}{1}(sel);
    norms = collection{i}{2}(sel);
    c_parameters = collection{i}{3}(sel);
    stability = zeros(1, length(sel));
    
    params = [omega T(i) P0 sigma1];
    for j = 1:length(sel)
        fprintf('--> %g/%g\n', j, length(sel))
        params(1) = omegas(j);
        [u0, du0] = asympt_left_cosine_nho(params, c_parameters(j), xspan(1));
        [Xh, Uh] = f_solve_cosine_nho(params, xspan, [u0 du0], 8192);
        if strcmp(type, 'even')
            X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); Uh((end-1):-1:1, 1)];
        else
            X = [Xh; -Xh((end-1):-1:1)]; U = [Uh(:, 1); -Uh((end-1):-1:1, 1)];
        end
        [stability(j), max_real, min_abs] = is_stable_nho(params, X, U, 128);
        fprintf('max: %g, min: %g\n', max_real, min_abs)
    end
    
    reduced{i} = {omegas, norms, c_parameters, stability};
end

reduced_G3 = reduced;

figure('Position', [100, 100, 400, 200]); hold on
axis([-1 8 0 60])
xticks(-1:1:8)

for i = 1:length(reduced_G0)
    plot_with_stability_smart(reduced_G0{i}{1}, reduced_G0{i}{2}, reduced_G0{i}{4});
end

for i = 1:length(reduced_G1)
    plot_with_stability_smart(reduced_G1{i}{1}, reduced_G1{i}{2}, reduced_G1{i}{4});
end

for i = 1:length(reduced_G2)
    plot_with_stability_smart(reduced_G2{i}{1}, reduced_G2{i}{2}, reduced_G2{i}{4});
end

for i = 1:length(reduced_G3)
    plot_with_stability_smart(reduced_G3{i}{1}, reduced_G3{i}{2}, reduced_G3{i}{4});
end

figure('Position', [100, 100, 400, 200]); hold on
axis([-1 8 0 6])
xticks(-2:1:8)

for i = 1:length(reduced_G0)
    plot_with_stability_smart(reduced_G0{i}{1}, i * ones(1, length(reduced_G0{i}{2})), reduced_G0{i}{4});
end

for i = 1:length(reduced_G1)
    sel = reduced_G1{i}{1} > 1.5;
    plot_with_stability_smart(reduced_G1{i}{1}(sel), i * ones(1, length(reduced_G1{i}{2}(sel))), reduced_G1{i}{4}(sel));
end

for i = 1:length(reduced_G2)
    sel = reduced_G2{i}{1} > 3.5;
    plot_with_stability_smart(reduced_G2{i}{1}(sel), i * ones(1, length(reduced_G2{i}{2}(sel))), reduced_G2{i}{4}(sel));
end

for i = 1:length(reduced_G3)
    sel = reduced_G3{i}{1} > 5.5;
    plot_with_stability_smart(reduced_G3{i}{1}(sel), i * ones(1, length(reduced_G3{i}{2}(sel))), reduced_G3{i}{4}(sel));
end