function  u_end = get_u_end_cosine(params, xspan, c, varargin)

if length(varargin) >= 1
	intervals = varargin{1};
else
    intervals = 8192;
end

u0 = [c*exp(sqrt(params(1)) * xspan(1)); c*exp(sqrt(params(1)) * xspan(1))];
[~,U] = f_solve_cosine(params, xspan, u0, intervals);
u_end = U(end, 1);

end

