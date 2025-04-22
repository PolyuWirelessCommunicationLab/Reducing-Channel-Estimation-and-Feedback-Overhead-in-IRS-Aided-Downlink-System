function vec = nearest_integers_avg(B, L)

    if B < 1 || L < 1 || floor(L) ~= L
        error('B must be a real number greater than 1 and L must be an integer greater than 1.');
    end
    
    tolerance = 0.01;
    max_iter = 100;
    
    lower_int = floor(B);    
    upper_int = ceil(B);
    
    if B == lower_int
        vec = repmat(B, 1, L);
        return;
    end
    
    vec = zeros(1, L);

    ratio = (B - lower_int) / (upper_int - lower_int);
    num_upper = round(ratio * L);
    num_lower = L - num_upper;
    
    avg = (num_upper * upper_int + num_lower * lower_int) / L;
    iter = 0;
    while abs(B - avg)/B > tolerance && iter <= max_iter
        if avg < B
            num_lower = num_lower - 1;
            num_upper = num_upper + 1;
        else
            num_lower = num_lower + 1;
            num_upper = num_upper - 1;
        end
        avg = (num_upper * upper_int + num_lower * lower_int) / L;
        iter = iter + 1;
    end
    
    vec(1:num_upper) = upper_int;
    vec(num_upper + 1:end) = lower_int;
    
    vec = vec(randperm(L));
end