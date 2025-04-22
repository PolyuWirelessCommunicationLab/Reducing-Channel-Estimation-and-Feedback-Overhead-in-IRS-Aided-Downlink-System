function distances = generate_user_positions(N)
    % Calculate the coordinates of the center of the circle
    x_c = (105^2 + 100^2 - 10^2) / 200;
    y_c_squared = 105^2 - x_c^2;
    y_c_positive = sqrt(y_c_squared);
    y_c_negative = -sqrt(y_c_squared);

    % select positive y_c or negative one
    circle_center_X = x_c;
    circle_center_Y = y_c_positive; % or y_c_negative
    
    % radius of the area of users
    circle_radius = 5;
    
    distances = zeros(1, N);
    
    for i = 1:N
        % generate user positions
        angle = 2 * pi * rand();
        r = circle_radius * sqrt(rand());
        user_X = circle_center_X + r * cos(angle);
        user_Y = circle_center_Y + r * sin(angle);
        
        % calculate distances between users and the IRS
        distances(i) = sqrt((user_X - 100)^2 + (user_Y - 0)^2);
    end
end
