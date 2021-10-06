function pairs = get_dichotomy_pairs(parameters, values)

pairs = {};

prev = sign(values(1));
for i = 2:length(values)
    curr = sign(values(i));
    
    if (prev == +1 && curr == -1) || (prev == -1 && curr == +1)
        pairs{end + 1} = [parameters(i - 1) parameters(i)];
    end
        
    prev = curr;
end

end
