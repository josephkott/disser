function [u, du] = asympt_left_cosine_nho(params, c, x)
% Asymtotical behaviour on the -\infty
%
% INPUT:
%	c - constant of the asymptote
%   x < 0

omega = params(1);

u  = c * (-x)^(0.5 * (omega - 1)) * exp(-(x^2) / 2);
du = -0.5 * c * (omega - 1) * exp(-(x^2) / 2) * (-x)^(0.5 * (omega - 1) - 1) - ...
    c * exp(-(x^2) / 2) * x * ((-x)^(0.5 * (omega - 1)));

end