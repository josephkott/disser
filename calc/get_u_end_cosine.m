function u_end = get_u_end_cosine(params, xspan, c, varargin)

if length(varargin) >= 1
	intervals = varargin{1};
else
    intervals = 8192;
end

init = [c*exp(sqrt(params(1)) * xspan(1)); c*exp(sqrt(params(1)) * xspan(1))];
[~, U] = f_solve_cosine(params, xspan, init, intervals);
u_end = U(end, 1);

end
