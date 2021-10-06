function ux_end = get_ux_end_cosine_nho(params, xspan, c, varargin)

if length(varargin) >= 1
	intervals = varargin{1};
else
    intervals = 8192;
end

[u0, du0] = asympt_left_cosine_nho(params, c, xspan(1));
[~, U] = f_solve_cosine_nho(params, xspan, [u0 du0], intervals);
ux_end = U(end, 2);

end

