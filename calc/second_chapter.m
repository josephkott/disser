%% Picture (1). Island set emergence
clc; clear

addpath('./mex/');

params = [1 0];
Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Up, true, 'gray')

%% Picture (1). Island set emergence, (a) no

M = 1024; N = 1024;
params = [1 0.6]; u_span = [-8 8]; du_span = [-40 40];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
figure('Position', [100, 100, 350, 220])
plot_scan(u_span, du_span, Up, false, 'gray')

%% Picture (1). Island set emergence, (b) no

M = 1024; N = 1024;
params = [1 0.3]; u_span = [-6 6]; du_span = [-18 18];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
figure('Position', [100, 100, 350, 220])
plot_scan(u_span, du_span, Up, false, 'gray')

%% Picture (1). Island set emergence, (c) no

M = 1024; N = 1024;
params = [1 0.1]; u_span = [-7 7]; du_span = [-20 20];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
figure('Position', [100, 100, 350, 220])
plot_scan(u_span, du_span, Up, false, 'gray')

%% Picture (1). Island set emergence, (d) yes

M = 1024; N = 1024;
params = [1 -0.1]; u_span = [-8.5 8.5]; du_span = [-25 25];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
figure('Position', [100, 100, 350, 220])
plot_scan(u_span, du_span, Up, false, 'gray')

%% Picture (2). Hypotheses for the island set, $\omega = 1.5$, $\alpha = 0$
clc; clear

addpath('./mex/')
params = [1.5 0];

%% Picture (2). Main plot

M = 2048; N = 2048;
u_span = [-8 8]; du_span = [-22 22];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
figure('Position', [100, 100, 350, 220])
plot_scan(u_span, du_span, Up, false, 'gray')
xticks(-8:8)
yticks(-20:5:20)

%% Picure (2). Island "0", "+" set

M = 1024; N = 1024;
u_span = [-0.8 0.8]; du_span = [-1 1];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Up, true, 'gray')

%% Picture (2). Island "0", "-" set

M = 1024; N = 1024;
u_span = [-0.8 0.8]; du_span = [-1 1];

Um = f_scan_cosine(params, [0, -pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Um, true, 'gray')

%% Picture (2). Island "0", "2+" set

M = 1024; N = 1024;
u_span = [-0.8 0.8]; du_span = [-1 1];

Up = f_scan_bounded_cosine(params, [0, pi], u_span, du_span, M, N, 4, 20, 1024);
U2p = f_scan_cosine(params, [0, 2*pi], u_span, du_span, M, N, 1024);

U2pb = ones(N, M);
U2pb(Up == 0 & U2p == 0) = 0;

plot_scan(u_span, du_span, U2pb, true, 'gray')

%% Picture (2). Island "0", "2-" set

M = 1024; N = 1024;
u_span = [-0.8 0.8]; du_span = [-1 1];

Um = f_scan_bounded_cosine(params, [0, -pi], u_span, du_span, M, N, 4, 20, 1024);
U2m = f_scan_cosine(params, [0, -2*pi], u_span, du_span, M, N, 1024);

U2mb = ones(N, M);
U2mb(Um == 0 & U2m == 0) = 0;

plot_scan(u_span, du_span, U2mb, true, 'gray')

%% Picture (2). Island "0", signs of the $D \mathcal{P}^{-1}$ operator (H)

M = 128; N = 128;
u_span = [-0.8 0.8]; du_span = [-1 1];

backward_linearization_plots_cosine(params, u_span, du_span, M, N)

%% Picture (2). Island "0", signs of the $D \mathcal{P}$ operator (V)

M = 128; N = 128;
u_span = [-0.8 0.8]; du_span = [-1 1];

forward_linearization_plots_cosine(params, u_span, du_span, M, N)

%% Picure (2). Island "+1", "+" set

M = 1024; N = 1024;
u_span = [1.6 2.8]; du_span = [-3 3];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Up, true, 'gray')

%% Picture (2). Island "+1", "-" set

M = 1024; N = 1024;
u_span = [1.6 2.8]; du_span = [-3 3];

Um = f_scan_cosine(params, [0, -pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Um, true, 'gray')

%% Picture (2). Island "+1", "2+" set

M = 1024; N = 1024;
u_span = [1.6 2.8]; du_span = [-3 3];

Up = f_scan_bounded_cosine(params, [0, pi], u_span, du_span, M, N, 4, 20, 1024);
U2p = f_scan_cosine(params, [0, 2*pi], u_span, du_span, M, N, 1024);

U2pb = ones(N, M);
U2pb(Up == 0 & U2p == 0) = 0;

plot_scan(u_span, du_span, U2pb, true, 'gray')

%% Picture (2). Island "+1", "2-" set

M = 1024; N = 1024;
u_span = [1.6 2.8]; du_span = [-3 3];

Um = f_scan_bounded_cosine(params, [0, -pi], u_span, du_span, M, N, 4, 20, 1024);
U2m = f_scan_cosine(params, [0, -2*pi], u_span, du_span, M, N, 1024);

U2mb = ones(N, M);
U2mb(Um == 0 & U2m == 0) = 0;

plot_scan(u_span, du_span, U2mb, true, 'gray')

%% Picture (2). Island "+1", signs of the $D \mathcal{P}^{-1}$ operator (H)

M = 1024; N = 1024;
u_span = [1.6 2.8]; du_span = [-3 3];

backward_linearization_plots_cosine(params, u_span, du_span, M, N)

% From the up to down: 7 2 8 2 7

%% Picture (2). Island "+1", signs of the $D \mathcal{P}$ operator (V)

M = 1024; N = 1024;
u_span = [1.6 2.8]; du_span = [-3 3];

forward_linearization_plots_cosine(params, u_span, du_span, M, N)

% From the up to down: 3 6 4 6 3

%% Picure (2). Island "+1i", "+" set

M = 1024; N = 1024;
u_span = [-0.4 +0.4]; du_span = [12.5 14.3];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Up, true, 'gray')

%% Picture (2). Island "+1i", "-" set

M = 1024; N = 1024;
u_span = [-0.4 +0.4]; du_span = [12.5 14.3];

Um = f_scan_cosine(params, [0, -pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Um, true, 'gray')

%% Picture (2). Island "+1i", "2+" set

M = 1024; N = 1024;
u_span = [-0.4 +0.4]; du_span = [12.5 14.3];

Up = f_scan_bounded_cosine(params, [0, pi], u_span, du_span, M, N, 4, 20, 1024);
U2p = f_scan_cosine(params, [0, 2*pi], u_span, du_span, M, N, 1024);

U2pb = ones(N, M);
U2pb(Up == 0 & U2p == 0) = 0;

plot_scan(u_span, du_span, U2pb, true, 'gray')

%% Picture (2). Island "+1i", "2-" set

M = 1024; N = 1024;
u_span = [-0.4 +0.4]; du_span = [12.5 14.3];

Um = f_scan_bounded_cosine(params, [0, -pi], u_span, du_span, M, N, 4, 20, 1024);
U2m = f_scan_cosine(params, [0, -2*pi], u_span, du_span, M, N, 1024);

U2mb = ones(N, M);
U2mb(Um == 0 & U2m == 0) = 0;

plot_scan(u_span, du_span, U2mb, true, 'gray')

%% Picture (2). Island "+1i", signs of the $D \mathcal{P}^{-1}$ operator (H)

M = 1024; N = 1024;
u_span = [-0.4 +0.4]; du_span = [12.5 14.3];

backward_linearization_plots_cosine(params, u_span, du_span, M, N)

%% Picture (2). Island "+1i", signs of the $D \mathcal{P}$ operator (V)

M = 1024; N = 1024;
u_span = [-0.4 +0.4]; du_span = [12.5 14.3];

forward_linearization_plots_cosine(params, u_span, du_span, M, N)

%% Picure (2). Island "+2", "+" set

M = 1024; N = 1024;
u_span = [7.06 7.16]; du_span = [-1.2 1.2];

Up = f_scan_cosine(params, [0, pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Up, true, 'gray')

%% Picture (2). Island "+2", "-" set

M = 1024; N = 1024;
u_span = [7.06 7.16]; du_span = [-1.2 1.2];

Um = f_scan_cosine(params, [0, -pi], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Um, true, 'gray')

%% Picture (2). Island "+2", "2+" set

M = 1024; N = 1024;
u_span = [7.06 7.16]; du_span = [-1.2 1.2];

Up = f_scan_bounded_cosine(params, [0, pi], u_span, du_span, M, N, 4, 20, 1024);
U2p = f_scan_cosine(params, [0, 2*pi], u_span, du_span, M, N, 1024);

U2pb = ones(N, M);
U2pb(Up == 0 & U2p == 0) = 0;

plot_scan(u_span, du_span, U2pb, true, 'gray')

%% Picture (2). Island "+2", "2-" set

M = 1024; N = 1024;
u_span = [7.06 7.16]; du_span = [-1.2 1.2];

Um = f_scan_bounded_cosine(params, [0, -pi], u_span, du_span, M, N, 4, 20, 1024);
U2m = f_scan_cosine(params, [0, -2*pi], u_span, du_span, M, N, 1024);

U2mb = ones(N, M);
U2mb(Um == 0 & U2m == 0) = 0;

plot_scan(u_span, du_span, U2mb, true, 'gray')

%% Picture (2). Island "+2", signs of the $D \mathcal{P}^{-1}$ operator (H)

M = 1024; N = 1024;
u_span = [7.06 7.16]; du_span = [-1.2 1.2];

backward_linearization_plots_cosine(params, u_span, du_span, M, N)

%% Picture (2). Island "+2", signs of the $D \mathcal{P}$ operator (V)

M = 1024; N = 1024;
u_span = [7.06 7.16]; du_span = [-1.2 1.2];

forward_linearization_plots_cosine(params, u_span, du_span, M, N)

%% Picture (3). Solutions for cosine equation, pi-periodic solution
% Code: (+1 +1 +1 +1 +1 +1)
clc; clear
addpath('./mex/')

params = [1.5 0];
L = pi; xspan = [0 pi];

f_map = @(p) f_map_cosine(params, p, 2^16);
problem = @(p) p - f_map(p);
solution = fsolve(problem, [2, -0.45]);
[X, U] = f_solve_cosine(params, xspan, solution, 2^16);

figure('Position', [100, 100, 350, 200])
X_plot = [-X(end:-1:1); X; X + L; X + 2*L; X + 3*L; X + 4*L; X + 5*L];
U_plot = [U(:, 1); U(:, 1); U(:, 1); U(:, 1); U(:, 1); U(:, 1); U(:, 1)];
plot(X_plot, U_plot, 'Color', 'black')

axis([-L 6*L 0 2.5])
xticks([0 L 2*L 3*L 4*L 5*L])
yticks([0 1 2])
xticklabels({'0','+\pi','+2\pi','+3\pi', '+4\pi', '+5\pi'})
grid on

%% Picture (3). Solutions for cosine equation, 2pi-periodic solution
% Code: (+1 -1 +1 -1 +1 -1)
clc; clear
addpath('./mex/')

params = [1.5 0];
L = pi; xspan = [0 2*pi];

f_map = @(p) f_map_cosine(params, p, 2^18);
problem = @(p) p - f_map(f_map(p));
solution = fsolve(problem, [2.25103 0.120235]);
[X, U] = f_solve_cosine(params, xspan, solution, 2*2^18);

figure('Position', [100, 100, 350, 200])
X_plot = [-X(end:-1:1); X; X + 2*L; X + 4*L];
U_plot = [U(:, 1); U(:, 1); U(:, 1); U(:, 1)];
plot(X_plot, U_plot, 'Color', 'black')

axis([-L 6*L -2.5 2.5])
xticks([0 L 2*L 3*L 4*L 5*L])
xticklabels({'0','+\pi','+2\pi','+3\pi', '+4\pi', '+5\pi'})
grid on

%% Picture (3). Solutions for cosine equation, another 2pi-periodic solution
% Code: (+1 0 +1 0 +1 0)

clc; clear
addpath('./mex/')

params = [1.5 0];
L = pi; xspan = [0 2*pi];

f_map = @(p) f_map_cosine(params, p, 2^18);
problem = @(p) p - f_map(f_map(p));
solution = fsolve(problem, [0.732747 -0.824047]);
[X, U] = f_solve_cosine(params, xspan, solution, 2*2^18);

figure('Position', [100, 100, 350, 200])
X_plot = [-X(end:-1:1); X; X + 2*L; X + 4*L];
U_plot = [U(:, 1); U(:, 1); U(:, 1); U(:, 1)];
plot(X_plot, U_plot, 'Color', 'black')

axis([-L 6*L 0 2.5])
xticks([0 L 2*L 3*L 4*L 5*L])
xticklabels({'0','+\pi','+2\pi','+3\pi', '+4\pi', '+5\pi'})
grid on

%% Picture (3). Solutions for cosine equation, localized FS
% Code: (0 0 +1 0 0)

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_ux_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 2;
c_right = 5;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
plot([X, -X(end:-1:1)], [U(:, 1) U(end:-1:1, 1)], 'Color', 'black')

axis([-3*pi 3*pi -0.5 2.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

%% Picture (3). Solutions for cosine equation, localized DS
% Code: (0 0 +1i 0 0)

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_u_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 4;
c_right = 8;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
plot([X, -X(end:-1:1)], [U(:, 1) -U(end:-1:1, 1)], 'Color', 'black')

axis([-3*pi 3*pi -5.5 5.5])
xticks([-2*pi -pi 0 pi 2*pi])
yticks([-4 -2 0 2 4])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

%% Picture (3). Solutions for cosine equation, (0 0 +2 0 0)

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_ux_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 9;
c_right = 10;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
X = [X; -X((end-1):-1:1)];
U = [-U(:, 1); -U((end-1):-1:1, 1)];
plot(X, U, 'Color', 'black')

axis([-3*pi 3*pi -9 9])
xticks([-2*pi -pi 0 pi 2*pi])
yticks([-8 -4 0 4 8])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

%% Picture (3). Solutions for cosine equation, (0 +1 0 +1 0)

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_ux_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 136.62; c_right = 138.45; % (... 0 +1 0 +1 0 ...)
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
X = [X; -X((end-1):-1:1)];
U = [U(:, 1); U((end-1):-1:1, 1)];
plot(X, U, 'Color', 'black')

axis([-3*pi 3*pi -0.5 2.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

%% Picture (3). Solutions for cosine equation, (0 +1 0 -1 0)

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_u_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 136.38; c_right = 144.35;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
plot([X, -X(end:-1:1)], [U(:, 1) -U(end:-1:1, 1)], 'Color', 'black')

axis([-3*pi 3*pi -2.5 2.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

%% Picture (3). Solutions for cosine equation, (0 +1 -1i +1i +1 0)

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-6 * pi,  pi / 2];

eps = 1e-9;
f = @(c) get_ux_end_cosine(params, xspan, c);
c_left = 128.5; c_right = 128.99;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
X = [X; -X(end:-1:1) + pi];
U = [U(:, 1); U(end:-1:1, 1)];
plot(X, U, 'Color', 'black')

axis([-3*pi 4*pi -5.5 5.5])
xticks([-2*pi -pi 0 pi 2*pi 3*pi])
yticks([-4 -2 0 2 4])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi', '3\pi'})
grid on

%% Picture (4) Stability demonstration, left column, (0 +1 0 -1 0) stable

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_u_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 136.38; c_right = 144.35;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
X = [X; -X((end-1):-1:1)];
U = [U(:, 1); -U((end-1):-1:1, 1)];
plot(X, U, 'Color', 'black')

axis([-3*pi 3*pi -2.5 2.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

eigenvalues = get_spectrum(params, X, U, 256);
figure('Position', [100, 100, 350, 250])
plot_spectrum(params, eigenvalues);
axis([-0.3 0.3 -3 3])

[Grid, Solution, Norm] = CFDS(params, X, U, 100, 0.05);

figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-3*pi 3*pi])
yticks([-2*pi -pi 0 pi 2*pi])
yticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})

%% Picture (4) Stability demonstration, middle column, (0 +1 +1 +1 0) unstable

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_ux_end_cosine(params, xspan, c);

eps = 1e-9;
% c_left = 139.15; c_right = 155.62; % (0 +1 -1 +1 0)
c_left = 126.92; c_right = 135.68; % (0 +1 +1 +1 0)
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
X = [X; -X((end-1):-1:1)];
U = [U(:, 1); U((end-1):-1:1, 1)];
plot(X, U, 'Color', 'black')

axis([-3*pi 3*pi -0.5 2.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

eigenvalues = get_spectrum(params, X, U, 512);
figure('Position', [100, 100, 350, 250])
plot_spectrum(params, eigenvalues);
axis([-1.5 1.5 -5 5])
yticks([-4 -2 0 2 4])
  
[Grid, Solution, Norm] = CFDS(params, X, U, 20, 0.02);

figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-3*pi 3*pi])
yticks([-2*pi -pi 0 pi 2*pi])
yticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})

%% Picture (4) Stability demonstration, middle column, (0 0 +2 0 0) unstable

clc; clear
addpath('./mex/')

params = [1.5 0];
xspan = [-5 * pi,  0];

f = @(c) get_ux_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 9;
c_right = 10;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 200])
X = [X; -X((end-1):-1:1)];
U = [-U(:, 1); -U((end-1):-1:1, 1)];
plot(X, U, 'Color', 'black')

axis([-3*pi 3*pi -9 9])
xticks([-2*pi -pi 0 pi 2*pi])
yticks([-8 -4 0 4 8])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

eigenvalues = get_spectrum(params, X, U, 512);
figure('Position', [100, 100, 350, 250])
plot_spectrum(params, eigenvalues);
axis([-3 3 -42 42])

[Grid, Solution, Norm] = CFDS(params, X, U, 10, 0.005);

figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-2*pi 2*pi])
yticks([-pi 0 pi])
yticklabels({'-\pi','0','\pi'})
