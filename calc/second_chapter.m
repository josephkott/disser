%% TODO: Copy all other stuff to the repository
% ...

%% Illustration to the Numerical Test of Existence of Coding Homeomorphism
addpath('./mex/')

clc; clear

L1 = 2; L2 = 1;
L = L1 + L2;
params = [1 L1 L2];

u_span = [-2 +2]; du_span = [-4 +4];
M = 1024; N = 1024;

Up = f_scan_piecewise(params, [0, L1], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Up, true, 'gray')

Um = f_scan_piecewise(params, [0 -L], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, Um, true, 'gray')

Up = f_scan_bounded_piecewise(params, [0, L1], u_span, du_span, M, N, 2, 4, 1024);
U2p = f_scan_piecewise(params, [0, L + L1], u_span, du_span, M, N, 1024);
U2pb = ones(N, M);
U2pb(Up == 0 & U2p == 0) = 0;
plot_scan(u_span, du_span, U2pb, true, 'gray')

Um = f_scan_bounded_piecewise(params, [0, -L], u_span, du_span, M, N, 2, 4, 1024);
U2m = f_scan_piecewise(params, [0, -2*L], u_span, du_span, M, N, 1024);
U2mb = ones(N, M);
U2mb(Um == 0 & U2m == 0) = 0;
plot_scan(u_span, du_span, U2mb, true, 'gray')

backward_linearization_plots_piecewise(params, u_span, du_span, M, N) % $D \mathcal{P}^{-1}$, (H)
forward_linearization_plots_piecewise(params, u_span, du_span, M, N) % $D \mathcal{P}$, (V)

%% Picture (2). Island "0", signs of the $D \mathcal{P}$ operator (V)

M = 128; N = 128;
u_span = [-0.8 0.8]; du_span = [-1 1];

forward_linearization_plots_cosine(params, u_span, du_span, M, N)