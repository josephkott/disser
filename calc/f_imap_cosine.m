function q = f_imap_cosine(params, p, intervals)

[~, U, is_inf] = f_solve_cosine(params, [0 -pi], p, intervals);

if is_inf
    q = NaN;
else
    q = U(end, :);
end

end

