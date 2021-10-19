function type = get_signs_configuration_type(J)

equals = @(x, y) all(all(x == y));

S = sign(J);
if equals(S, [+1 +1; +1 +1])
    type = 1;
elseif equals(S, [-1 -1; -1 -1])
    type = 2;
elseif equals(S, [+1 +1; -1 -1])
    type = 3;
elseif equals(S, [-1 -1; +1 +1])
    type = 4;
elseif equals(S, [+1 -1; -1 +1])
    type = 5;
elseif equals(S, [-1 +1; +1 -1])
    type = 6;
elseif equals(S, [+1 -1; +1 -1])
    type = 7;
elseif equals(S, [-1 +1; -1 +1])
    type = 8;
else
    type = -1;
    
end
