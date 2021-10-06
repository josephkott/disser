function u_end = get_u_end_cosine_nho(params, xspan, c, varargin)

if length(varargin) >= 1
	intervals = varargin{1};
else
    intervals = 8192;
end

[u0, du0] = asympt_left_cosine_nho(params, c, xspan(1));
[~, U] = f_solve_cosine_nho(params, xspan, [u0 du0], intervals); 
u_end = U(end, 1);

end

