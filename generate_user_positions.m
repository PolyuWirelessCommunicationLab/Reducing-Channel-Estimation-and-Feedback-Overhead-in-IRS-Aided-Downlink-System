function distances = generate_user_positions(N)
    % 计算圆心坐标
    x_c = (105^2 + 100^2 - 10^2) / 200;
    y_c_squared = 105^2 - x_c^2;
    y_c_positive = sqrt(y_c_squared); % 正的y_c值
    y_c_negative = -sqrt(y_c_squared); % 负的y_c值

    % 选择正的y_c值或者负的y_c值
    circle_center_X = x_c;
    circle_center_Y = y_c_positive; % 或者 y_c_negative，取决于你的需求
    
    % 用户区域的半径
    circle_radius = 5;
    
    % 初始化距离数组
    distances = zeros(1, N);
    
    for i = 1:N
        % 生成用户位置
        angle = 2 * pi * rand();
        r = circle_radius * sqrt(rand());
        user_X = circle_center_X + r * cos(angle);
        user_Y = circle_center_Y + r * sin(angle);
        
        % 计算用户与IRS之间的距离
        distances(i) = sqrt((user_X - 100)^2 + (user_Y - 0)^2);
    end
end