function  ux_end = get_ux_end_cosine(params, xspan, c)

u0 = [c*exp(sqrt(params(1)) * xspan(1)); c*exp(sqrt(params(1)) * xspan(1))];
[~,U] = f_solve_cosine(params, xspan, u0, 8192);
ux_end = U(end, 2);

end

