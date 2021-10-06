function [ stable, max_real, min_abs ] = is_stable_nho( params, X, U, N )
% Use @get_spectrum_nho to determine the stability of the mode
%
% INPUT:
% :params: parameter of NLS with nonlinear potential (see info.txt)
% :X: grid
% :U: values on grid
% :N: number of harmonics in fourie series

% OUTPUT:
% :stable: true of false
%

eigenvalues = get_spectrum_nho(params, X, U, N);
max_real = max(abs(real(eigenvalues)));
min_abs = min(abs(eigenvalues));

% Compareson with zero
if max_real < 1e-3
	stable = true;
else
	stable = false;
end

end
