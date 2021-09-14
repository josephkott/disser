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

%% Picture (4). Stability demonstration, left column, (0 +1 0 -1 0) stable

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

%% Picture (4). Stability demonstration, middle column, (0 +1 +1 +1 0) unstable

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

%% Picture (4). Stability demonstration, middle column, (0 0 +2 0 0) unstable

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

%% Picture (5). Branches with stability

clc; clear
addpath('./mex/')

params = [6 0];
xspan = [-5 * pi 0];

%% Picutre (5). Get FS

f = @(c) get_ux_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 5;
c_right = 10;
FS_c_root = dichotomy(f, c_left, c_right, eps);

u0 = FS_c_root * exp(sqrt(params(1)) * xspan(1));
du0 = FS_c_root * exp(sqrt(params(1)) * xspan(1));
[Xh_FS, Uh_FS] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 16384);

X_FS = [Xh_FS; -Xh_FS((end-1):-1:1)];
U_FS = [Uh_FS(:, 1); Uh_FS((end-1):-1:1, 1)];
figure; plot(X_FS, U_FS)

%% Picture (5). Continuation of FS

FS_omegas  = [6:-0.01:0.1, 0.099:-0.001:0.05];
FS_norms   = zeros(1, length(FS_omegas));
FS_c_roots = zeros(1, length(FS_omegas));
FS_stability = zeros(1, length(FS_omegas));

FS_c_roots(1) = FS_c_root;
FS_norms(1) = trapz(X_FS, U_FS .^ 2);
FS_stability(1) = is_stable([FS_omegas(1) 0], X_FS, U_FS, 128);

for i = 2:length(FS_omegas)
    f = @(c) get_ux_end_cosine([FS_omegas(i), 0], xspan, c);
    FS_c_roots(i) = newton(f, FS_c_roots(i - 1));
    
    u0 = FS_c_roots(i) * exp(sqrt(FS_omegas(i)) * xspan(1));
    du0 = FS_c_roots(i) * exp(sqrt(FS_omegas(i)) * xspan(1));
    [Xh_FS, Uh_FS] = f_solve_cosine([FS_omegas(i), 0], xspan, [u0 du0], 8192);
    
    X_FS = [Xh_FS; -Xh_FS((end-1):-1:1)];
    U_FS = [Uh_FS(:, 1); Uh_FS((end-1):-1:1, 1)];
    FS_norms(i) = trapz(X_FS, U_FS .^ 2);
    FS_stability(i) = is_stable([FS_omegas(i) 0], X_FS, U_FS, 128);
    
    fprintf('iter: %g/%g, omega: %g, c_root: %g, norm: %g, is_stable: %1.g\n', ...
            i, length(FS_omegas), FS_omegas(i), FS_c_roots(i), FS_norms(i), FS_stability(i))
end

%% Picture (5). Get DS

f = @(c) get_u_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 25;
c_right = 30;
DS_c_root = dichotomy(f, c_left, c_right, eps);

u0 = DS_c_root * exp(sqrt(params(1)) * xspan(1));
du0 = DS_c_root * exp(sqrt(params(1)) * xspan(1));
[Xh_DS, Uh_DS] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 16384);

X_DS = [Xh_DS; -Xh_DS((end-1):-1:1)];
U_DS = [Uh_DS(:, 1); -Uh_DS((end-1):-1:1, 1)];
figure; plot(X_DS, U_DS)

%% Picture (5). Continuation of DS

DS_omegas  = [6:-0.01:0.27, 0.2699:-0.0001:0.26];
DS_norms   = zeros(1, length(DS_omegas));
DS_c_roots = zeros(1, length(DS_omegas));
DS_stability = zeros(1, length(DS_omegas));
FS_stability(1) = is_stable([DS_omegas(1) 0], X_DS, U_DS, 256);

DS_c_roots(1) = DS_c_root;
DS_norms(1) = trapz(X_DS, U_DS .^ 2);

for i = 2:length(DS_omegas)
    f = @(c) get_u_end_cosine([DS_omegas(i), 0], xspan, c);
    DS_c_roots(i) = newton(f, DS_c_roots(i - 1));
    
    u0 = DS_c_roots(i) * exp(sqrt(DS_omegas(i)) * xspan(1));
    du0 = DS_c_roots(i) * exp(sqrt(DS_omegas(i)) * xspan(1));
    [Xh_DS, Uh_DS] = f_solve_cosine([DS_omegas(i), 0], xspan, [u0 du0], 8192);
    
    X_DS = [Xh_DS; -Xh_DS((end-1):-1:1)];
    U_DS = [Uh_DS(:, 1); -Uh_DS((end-1):-1:1, 1)];
    DS_norms(i) = trapz(X_DS, U_DS .^ 2);
    DS_stability(i) = is_stable([DS_omegas(i) 0], X_DS, U_DS, 256);
    
    fprintf('iter: %g/%g, omega: %g, c_root: %g, norm: %g, is_stable: %1.g\n', ...
            i, length(DS_omegas), DS_omegas(i), DS_c_roots(i), DS_norms(i), DS_stability(i))
end

%% Picture (5). Get (... 0 +1 -1i -1 0 ...), DSb (bifurcate with DS)

f = @(c) get_u_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 18391;
c_right = 18616;
DSb_c_root = dichotomy(f, c_left, c_right, eps);

u0 = DSb_c_root * exp(sqrt(params(1)) * xspan(1));
du0 = DSb_c_root * exp(sqrt(params(1)) * xspan(1));
[Xh_DSb, Uh_DSb] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 16384);

X_DSb = [Xh_DSb; -Xh_DSb((end-1):-1:1)];
U_DSb = [Uh_DSb(:, 1); -Uh_DSb((end-1):-1:1, 1)];
figure
plot(X_DSb, U_DSb)

%% Picture (5). Continuation of DSb

DSb_omegas  = [6:-0.005:0.27, 0.2699:-0.0001:0.26];
DSb_norms   = zeros(1, length(DSb_omegas));
DSb_c_roots = zeros(1, length(DSb_omegas));
DSb_stability = zeros(1, length(DSb_omegas));

DSb_c_roots(1) = DSb_c_root;
DSb_norms(1) = trapz(X_DSb, U_DSb .^ 2);
DSb_stability(1) = is_stable([DSb_omegas(1) 0], X_DSb, U_DSb, 256);

for i = 2:length(DSb_omegas)
    f = @(c) get_u_end_cosine([DSb_omegas(i), 0], xspan, c);
    DSb_c_roots(i) = newton(f, DSb_c_roots(i - 1));
    
    u0 = DSb_c_roots(i) * exp(sqrt(DSb_omegas(i)) * xspan(1));
    du0 = DSb_c_roots(i) * exp(sqrt(DSb_omegas(i)) * xspan(1));
    [Xh_DSb, Uh_DSb] = f_solve_cosine([DSb_omegas(i), 0], xspan, [u0 du0], 8192);
    
    X_DSb = [Xh_DSb; -Xh_DSb((end-1):-1:1)];
    U_DSb = [Uh_DSb(:, 1); -Uh_DSb((end-1):-1:1, 1)];
    DSb_norms(i) = trapz(X_DSb, U_DSb .^ 2);
    DSb_stability(i) = is_stable([DSb_omegas(i) 0], X_DSb, U_DSb, 256);
    
    fprintf('iter: %g/%g, omega: %g, c_root: %g, norm: %g, is_stable: %1.g\n', ...
            i, length(DSb_omegas), DSb_omegas(i), DSb_c_roots(i), DSb_norms(i), DSb_stability(i))
end

%% Picture (5). Get (0 +1 -1 +1 0), complex of FS

f = @(c) get_ux_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 18829;
c_right = 19044;
FSc_c_root = dichotomy(f, c_left, c_right, eps);

u0 = FSc_c_root * exp(sqrt(params(1)) * xspan(1));
du0 = FSc_c_root * exp(sqrt(params(1)) * xspan(1));
[Xh_FSc, Uh_FSc] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 16384);

X_FSc = [Xh_FSc; -Xh_FSc((end-1):-1:1)];
U_FSc = [Uh_FSc(:, 1); Uh_FSc((end-1):-1:1, 1)];
figure; plot(X_FSc, U_FSc)

%% Picture (5). Continuation of FSc

FSc_omegas  = [6:-0.005:0.1, 0.099:-0.001:0.05];
FSc_norms   = zeros(1, length(FSc_omegas));
FSc_c_roots = zeros(1, length(FSc_omegas));
FSc_stability = zeros(1, length(FSc_omegas));

FSc_c_roots(1) = FSc_c_root;
FSc_norms(1) = trapz(X_FSc, U_FSc .^ 2);
FSc_stability(1) = is_stable([FSc_omegas(1) 0], X_FSc, U_FSc, 256);

for i = 2:length(FSc_omegas)
    f = @(c) get_ux_end_cosine([FSc_omegas(i), 0], xspan, c);
    FSc_c_roots(i) = newton(f, FSc_c_roots(i - 1));
    
    u0 = FSc_c_roots(i) * exp(sqrt(FSc_omegas(i)) * xspan(1));
    du0 = FSc_c_roots(i) * exp(sqrt(FSc_omegas(i)) * xspan(1));
    [Xh_FSc, Uh_FSc] = f_solve_cosine([FSc_omegas(i), 0], xspan, [u0 du0], 8192);
    
    X_FSc = [Xh_FSc; -Xh_FSc((end-1):-1:1)];
    U_FSc = [Uh_FSc(:, 1); Uh_FSc((end-1):-1:1, 1)];
    FSc_norms(i) = trapz(X_FSc, U_FSc .^ 2);
    FSc_stability(i) = is_stable([FSc_omegas(i) 0], X_FSc, U_FSc, 256);
    
    fprintf('iter: %g/%g, omega: %g, c_root: %g, norm: %g, is_stable: %1.g\n', ...
            i, length(FSc_omegas), FSc_omegas(i), FSc_c_roots(i), FSc_norms(i), FSc_stability(i))
end

%% Picture (5). Get (0 +1 0 -1 0), 2nd complex of FS

f = @(c) get_u_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 18740;
c_right = 18877;
FSc2_c_root = dichotomy(f, c_left, c_right, eps);

u0 = FSc2_c_root * exp(sqrt(params(1)) * xspan(1));
du0 = FSc2_c_root * exp(sqrt(params(1)) * xspan(1));
[Xh_FSc2, Uh_FSc2] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 16384);

X_FSc2 = [Xh_FSc2; -Xh_FSc2((end-1):-1:1)];
U_FSc2 = [Uh_FSc2(:, 1); -Uh_FSc2((end-1):-1:1, 1)];
figure; plot(X_FSc2, U_FSc2)

%% Picture (5). Coninuation of FSc2

FSc2_omegas  = 6:-0.001:0.05;
FSc2_norms   = zeros(1, length(FSc2_omegas));
FSc2_c_roots = zeros(1, length(FSc2_omegas));
FSc2_stability = zeros(1, length(FSc2_omegas));

FSc2_c_roots(1) = FSc2_c_root;
FSc2_norms(1) = trapz(X_FSc2, U_FSc2 .^ 2);
FSc2_stability(1) = is_stable([FSc2_omegas(1) 0], X_FSc2, U_FSc2, 256);

for i = 2:length(FSc2_omegas)
    f = @(c) get_u_end_cosine([FSc2_omegas(i), 0], xspan, c);
    FSc2_c_roots(i) = newton(f, FSc2_c_roots(i - 1));
    
    u0 = FSc2_c_roots(i) * exp(sqrt(FSc2_omegas(i)) * xspan(1));
    du0 = FSc2_c_roots(i) * exp(sqrt(FSc2_omegas(i)) * xspan(1));
    [Xh_FSc2, Uh_FSc2] = f_solve_cosine([FSc2_omegas(i), 0], xspan, [u0 du0], 8192);
    
    X_FSc2 = [Xh_FSc2; -Xh_FSc2((end-1):-1:1)];
    U_FSc2 = [Uh_FSc2(:, 1); Uh_FSc2((end-1):-1:1, 1)];
    FSc2_norms(i) = trapz(X_FSc2, U_FSc2 .^ 2);
    FSc2_stability(i) = is_stable([FSc2_omegas(i) 0], X_FSc2, U_FSc2, 256);
    
    fprintf('iter: %g/%g, omega: %g, c_root: %g, norm: %g, is_stable: %1.g\n', ...
            i, length(FSc2_omegas), FSc2_omegas(i), FSc2_c_roots(i), FSc2_norms(i), FSc2_stability(i))
end

%% Picture (5). Wrap everything up

figure('Position', [100, 100, 350, 300])
hold on

plot_with_stability(FS_omegas, FS_norms, FS_stability)
plot_with_stability(DS_omegas(1:619), DS_norms(1:619), DS_stability(1:619))
plot_with_stability(DSb_omegas(1:1192), DSb_norms(1:1192), DSb_stability(1:1192))
plot_with_stability(FSc_omegas, FSc_norms, FSc_stability)
plot_with_stability(FSc2_omegas, FSc2_norms, FSc2_stability)

index = find(FS_omegas == 4);

u0 = FS_c_roots(index) * exp(sqrt(FS_omegas(index)) * xspan(1));
du0 = FS_c_roots(index) * exp(sqrt(FS_omegas(index)) * xspan(1));
[Xh_FS, Uh_FS] = f_solve_cosine([FS_omegas(index), 0], xspan, [u0 du0], 8192);
    
X_FS = [Xh_FS; -Xh_FS((end-1):-1:1)];
U_FS = [Uh_FS(:, 1); Uh_FS((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 250])
plot(X_FS, U_FS, 'Color', 'black')

axis([-3*pi 3*pi -0.5 3.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})

index = find(FSc2_omegas == 4);

u0 = FSc2_c_roots(index) * exp(sqrt(FSc2_omegas(index)) * xspan(1));
du0 = FSc2_c_roots(index) * exp(sqrt(FSc2_omegas(index)) * xspan(1));
[Xh_FSc2, Uh_FSc2] = f_solve_cosine([FSc2_omegas(index), 0], xspan, [u0 du0], 8192);
    
X_FSc2 = [Xh_FSc2; -Xh_FSc2((end-1):-1:1)];
U_FSc2 = [Uh_FSc2(:, 1); -Uh_FSc2((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 250])
plot(X_FSc2, U_FSc2, 'Color', 'black')

axis([-3*pi 3*pi -3.5 3.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})

index = find(FSc_omegas == 4);

u0 = FSc_c_roots(index) * exp(sqrt(FSc_omegas(index)) * xspan(1));
du0 = FSc_c_roots(index) * exp(sqrt(FSc_omegas(index)) * xspan(1));
[Xh_FSc, Uh_FSc] = f_solve_cosine([FSc_omegas(index), 0], xspan, [u0 du0], 8192);
    
X_FSc = [Xh_FSc; -Xh_FSc((end-1):-1:1)];
U_FSc = [Uh_FSc(:, 1); Uh_FSc((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 250])
plot(X_FSc, U_FSc, 'Color', 'black')

axis([-3*pi 3*pi -3.5 3.5])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})

index = find(DS_omegas == 4);

u0 = DS_c_roots(index) * exp(sqrt(DS_omegas(index)) * xspan(1));
du0 = DS_c_roots(index) * exp(sqrt(DS_omegas(index)) * xspan(1));
[Xh_DS, Uh_DS] = f_solve_cosine([DS_omegas(index), 0], xspan, [u0 du0], 8192);
    
X_DS = [Xh_DS; -Xh_DS((end-1):-1:1)];
U_DS = [Uh_DS(:, 1); -Uh_DS((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 250])
plot(X_DS, U_DS, 'Color', 'black')

axis([-3*pi 3*pi -6 6])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})

index = find(DSb_omegas == 4);

u0 = DSb_c_roots(index) * exp(sqrt(DSb_omegas(index)) * xspan(1));
du0 = DSb_c_roots(index) * exp(sqrt(DSb_omegas(index)) * xspan(1));
[Xh_DSb, Uh_DSb] = f_solve_cosine([DSb_omegas(index), 0], xspan, [u0 du0], 8192);
    
X_DSb = [Xh_DSb; -Xh_DSb((end-1):-1:1)];
U_DSb = [Uh_DSb(:, 1); -Uh_DSb((end-1):-1:1, 1)];

figure('Position', [100, 100, 350, 250])
plot(X_DSb, U_DSb, 'Color', 'black')

axis([-3*pi 3*pi -6 6])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})

%% Picture (6). Dipole soliton stability, stable FS

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

[Grid, Solution, Norm] = CFDS(params, X, U, 400, 0.05);

figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-2*pi 2*pi])
yticks([-pi 0 pi])
yticklabels({'-\pi','0','\pi'})

%% Picture (6). Dipole soliton stability, stable DS

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
X = [X; -X((end-1):-1:1)];
U = [U(:, 1); -U((end-1):-1:1, 1)];
plot(X, U, 'Color', 'black')

axis([-3*pi 3*pi -5.5 5.5])
xticks([-2*pi -pi 0 pi 2*pi])
yticks([-4 -2 0 2 4])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

eigenvalues = get_spectrum(params, X, U, 512);
figure('Position', [100, 100, 350, 250])
plot_spectrum(params, eigenvalues);
axis([-1.5 1.5 -5 5])
yticks([-4 -2 0 2 4])

[Grid, Solution, Norm] = CFDS(params, X, U, 400, 0.05);

figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-2*pi 2*pi])
yticks([-pi 0 pi])
yticklabels({'-\pi','0','\pi'})

%% Picture (6). Dipole soliton stability, unstable DS

clc; clear
addpath('./mex/')

params = [0.4 0];
xspan = [-5 * pi,  0];

f = @(c) get_u_end_cosine(params, xspan, c);

eps = 1e-9;
c_left = 2;
c_right = 4;
c_root = dichotomy(f, c_left, c_right, eps);

u0 = c_root * exp(sqrt(params(1)) * xspan(1));
du0 = c_root * exp(sqrt(params(1)) * xspan(1));
[X, U] = f_solve_cosine(params, [xspan(1) xspan(2)], [u0 du0], 8192);

X = [X; -X((end-1):-1:1)];
U = [U(:, 1); -U((end-1):-1:1, 1)];
figure('Position', [100, 100, 350, 200])
plot(X, U, 'Color', 'black')

axis([-3*pi 3*pi -5.5 5.5])
xticks([-2*pi -pi 0 pi 2*pi])
yticks([-4 -2 0 2 4])
xticklabels({'-2\pi','-\pi','0','\pi', '2\pi'})
grid on

eigenvalues = get_spectrum(params, X, U, 512);
figure('Position', [100, 100, 350, 250])
plot_spectrum(params, eigenvalues);
axis([-1.5 1.5 -5 5])
yticks([-4 -2 0 2 4])

[Grid, Solution, Norm] = CFDS(params, X, U, 400, 0.05);

figure('Position', [100, 100, 350, 250])
plot_evolution(Grid, Solution, [-2*pi 2*pi])
yticks([-pi 0 pi])
yticklabels({'-\pi','0','\pi'})

%% Usefull

c = 0:0.1:50;
left = zeros(length(c), 2);

for i = 1:length(c)
	u0 = c(i) * exp(sqrt(params(1)) * xspan(1));
	du0 = c(i) * exp(sqrt(params(1)) * xspan(1));
	[~, U] = f_solve_cosine(params, xspan, [u0 du0], 8192);
	left(i, :) = U(end, :);
end

plot(c, left(:, 1))