function x_root = newton( f, x0 )
%
% INPUT:
%	

MAX_ITERATION = 10; % 50

eps = 1e-9; % 1e-6
delta = 1e-6; % 1e-6

df = @(x) (f(x + eps) - f(x - eps)) / (2 * eps);

iter = 0;
while true
	iter = iter + 1;
	
    if iter > MAX_ITERATION
		x_root = NaN;
		return
    end
    
	x1 = x0 - f(x0) / df(x0);
    
	if abs(x1 - x0) < delta
        break
    else
        x0 = x1;
	end
end

x_root = x1;

end

