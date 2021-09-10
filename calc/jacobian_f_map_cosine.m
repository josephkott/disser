function J = jacobian_f_map_cosine(params, p, intervals)

eps = 1e-6;

p_left  = [p(1) - eps, p(2)];
p_right = [p(1) + eps, p(2)];
p_down  = [p(1), p(2) - eps];
p_up    = [p(1), p(2) + eps];

q_left  = f_map_cosine(params,  p_left, intervals);
q_right = f_map_cosine(params, p_right, intervals);
q_down  = f_map_cosine(params,  p_down, intervals);
q_up    = f_map_cosine(params,    p_up, intervals);

J = [(q_right(1) - q_left(1)) / (2 * eps), (q_right(2) - q_left(2)) / (2 * eps);
     (q_up(1) - q_down(1)) / (2 * eps), (q_up(2) - q_down(2)) / (2 * eps)];

end
