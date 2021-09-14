function  ux_end = get_ux_end_cosine(params, xspan, c, varargin)

if length(varargin) >= 1
	intervals = varargin{1};
else
    intervals = 8192;
end

u0 = [c*exp(sqrt(params(1)) * xspan(1)); c*exp(sqrt(params(1)) * xspan(1))];
[~,U] = f_solve_cosine(params, xspan, u0, intervals);
ux_end = U(end, 2);

end

