function q = f_map_piecewise(params, p, intervals)

L1 = params(2); L2 = params(3);
L = L1 + L2;

[~, U, is_inf] = f_solve_piecewise(params, [0 L], p, intervals);

if is_inf
    q = NaN;
else
    q = U(end, :);
end

end

