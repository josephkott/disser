function  u_end = get_u_end_cosine(params, xspan, c)

u0 = [c*exp(sqrt(params(1)) * xspan(1)); c*exp(sqrt(params(1)) * xspan(1))];
[~,U] = f_solve_cosine(params, xspan, u0, 8192);
u_end = U(end, 1);

end

