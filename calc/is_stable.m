function [ stable ] = is_stable( params, X, U, N )
% Use @get_spectrum to determine the stability of the mode
%
% INPUT:
% :params: parameter of NLS with nonlinear potential (see info.txt)
% :X: grid
% :U: values on grid
% :N: number of harmonics in fourie series

% OUTPUT:
% :stable: true of false
%

eigenvalues = get_spectrum(params, X, U, N);
max_real = max(abs(real(eigenvalues)));

% Compareson with zero
if max_real < 1e-2
	stable = true;
else
	stable = false;
end

% For debug purpoces
% plot_spectrum(params, eigenvalues);
% pause()
	
end

