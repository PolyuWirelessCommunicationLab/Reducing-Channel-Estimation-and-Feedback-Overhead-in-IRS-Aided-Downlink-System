function numbers = generate_vector(N, M)
    if M < N || M > 10 * N
        error('The sum M must be no less than N and no greater than 10 * N');
    end

    % Initialize the sequence with ones.
    numbers = ones(N,1);

    % Distribute the remaining part of M uniformly over the sequence.
    remaining = M - N;
    while remaining > 0
        % Randomly select an index in the sequence.
        index = randi([1, N]);

        % Calculate the possible increment such that the number does not exceed 12.
        increment = min(10 - numbers(index), remaining);

        % Randomly select an increment from [1, possible increment] and update the number.
        if increment > 1
            increment = randi([1, increment]);
        end
        numbers(index) = numbers(index) + increment;

        % Update the remaining part.
        remaining = remaining - increment;
    end
end