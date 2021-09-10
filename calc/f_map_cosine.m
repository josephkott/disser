function q = f_map_cosine(params, p, intervals)

if isnan(p)
    q = NaN;
    return
end

[~, U, is_inf] = f_solve_cosine(params, [0 pi], p, intervals);

if is_inf
    q = NaN;
else
    q = U(end, :);
end

end

