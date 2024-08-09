clear
clc

tic;

results = simulate_annealing_ris();
save('outage.mat', 'results');

toc;

function results = simulate_annealing_ris()
    % 参数设置
    % RIS coordinates
    % initial_ris_coords = [0, 0; 0, 15; 20, 0; 20, 15; 10, 5];
    initial_ris_coords = [12, 0; 0, 6; 15, 5];
    numRIS = size(initial_ris_coords, 1); % RIS数量
    Rx_trajectory = load('user_trajectories.mat');
    nunDevice = 6;
    total_elements = 128; % 每个RIS的元素数量
    max_iterations = 80; % 最大迭代次数
    T_init = 100; % 初始温度
    alpha = 0.99; % 降温速率
    S = [1, 10, 20, 30]; % 选取steps of users
    
    % 初始化存储结果的变量
    results = cell(length(S), 1);
    
    % 遍历每个时间片
    for t = 1:length(S)
        all_assignments = [];
        steps = S(t);
        Rx_coords = zeros(nunDevice, 2); % Rx coordinates
        for user =1:nunDevice
            Rx_coords(user,:) = Rx_trajectory.trajectories{user}(steps, :);
        end
        % 随机初始化设备和RIS的分配
        assignment = randi([1, numRIS], 1, nunDevice);
        [all_ris_coords, all_elements_per_ris, new_assignment] = update_elements_code(assignment, initial_ris_coords, total_elements);
        
        % 初始化中断概率
        [current_LoS, current_cost] = SA_Pout(Rx_coords, all_ris_coords, all_elements_per_ris, new_assignment);
        best_cost = current_cost;
        best_assignment = assignment;
        
        all_assignments = [all_assignments;assignment];

        T = T_init;
        cost_history = zeros(max_iterations, 1); % 初始化成本历史记录
        
        % 模拟退火过程
        for iter = 1:max_iterations
            % 生成新解
            up_assignment = generate_new_solution(assignment, numRIS);
            [all_ris_coords, all_elements_per_ris, new_assignment] = update_elements_code(up_assignment, initial_ris_coords, total_elements);
            
            % 判断新的assignment是否已经出现过            
            [isPresent, rowIndex] = ismember(up_assignment, all_assignments, 'rows');
            if isPresent
                % fprintf('Array is present in row %d\n', rowIndex);
                cost_history(iter) = best_cost;
                continue;
            else
                % fprintf('Array is not present in any row\n');
            end

            all_assignments = [all_assignments;up_assignment];
            [new_LoS, new_cost] = SA_Pout(Rx_coords, all_ris_coords, all_elements_per_ris, new_assignment);
            
            % 接受准则
            if new_cost < current_cost || exp((current_cost - new_cost) / T) > rand()
                assignment = up_assignment;
                current_cost = new_cost;
                
                if new_cost < best_cost
                    best_cost = new_cost;
                    best_assignment = up_assignment;
                end
            end
            
            % 记录当前最佳成本
            cost_history(iter) = best_cost;
            
            % 降温
            T = T * alpha;
            fprintf('Iter %d finished - Sum Outage Probability: %.4f\n', iter, best_cost);
        end
        
        % 存储结果
        results{t} = struct('assignment', best_assignment, 'cost', best_cost, 'cost_history', cost_history);
        fprintf('Time slot %d: Best cost = %.4f\n', t, best_cost);
    end
    
    % 绘制收敛图
    figure;
    axes2=axes('position',[0.2,0.1,0.74,0.7]); %这个是figure里面图的位置和大小，分别为离下边，左边的距离，还有图的高和宽
    set(gca, "LooseInset", [0,0,0,0]);%消除白边
    hold on;
    colors = {[193 018 033]/255, [120 000 001]/255, [102 155 187]/255, [000 047 073]/255};  % Define colors for different N values
    for t = 1:length(S)
        plot(results{t}.cost_history ./ 6, '-', 'LineWidth', 1.3, 'Color', colors{t}, 'DisplayName', ['Steps ', num2str(S(t))]);
    end
    % 
    % set(gca, 'YScale', 'log');
    % xlim([1, 80]);
    % ylim([0.001, 1]);
    grid on;
    set(gca, "FontSize",10, "Fontname", "Times new roman"); %轴刻度标签的字体大小和名称
    xlabel('Iteration',"FontSize",12, "Fontname", "Times new roman");
    ylabel('Outage Probability',"FontSize",12, "Fontname", "Times new roman");
    h = legend('show', "Fontname", "Times new roman", 'location','northeast', "FontSize",10, 'edgecolor', [1,1,1]);
    set(h, 'box', 'off'); %设置legend背景色透明
    hold off;
end

function [all_ris_coords, all_elements_per_ris, new_assignment] = update_elements_code(assignment, initial_ris_coords, total_elements)
% 存储每个小RIS的坐标和元素数量
    all_ris_coords = [];
    all_elements_per_ris = [];
    new_assignment = zeros(size(assignment));
    % 处理每个RIS编码
    unique_ris = unique(assignment);
    current_index = 1;
    for i = 1:length(unique_ris)
        ris_id = unique_ris(i);
        count = sum(assignment == ris_id);
    
        if count > 1
            % 将大RIS分割为更小的RIS
            [ris_coords, elements_per_ris] = splitRIS(initial_ris_coords(ris_id, :), total_elements, count);
        else
            % 单个RIS，不需要分割
            ris_coords = initial_ris_coords(ris_id, :);
            elements_per_ris = total_elements;
        end
        
        % 存储结果
        all_ris_coords = [all_ris_coords; ris_coords];
        all_ris_coords = abs(all_ris_coords);
        all_elements_per_ris = [all_elements_per_ris; repmat(elements_per_ris, count, 1)];
    
        % 更新device对应的RIS编码
        new_assignment(assignment == ris_id) = current_index:current_index + count - 1;
        current_index = current_index + count;
    end
end

function new_assignment = generate_new_solution(assignment, numRIS)
    % 随机选择一个设备并重新分配到另一个随机RIS
    new_assignment = assignment;
    n = length(assignment);
    % idx = randi([1, n]);
    % new_ris = randi([1, numRIS]);
    % new_assignment(idx) = new_ris;

    % 局部搜索，通过微调当前分配生成新解
    idx = randi([1, n]);
    neighbor_ris = mod(new_assignment(idx) + randi([-1, 1]), numRIS) + 1;
    new_assignment(idx) = neighbor_ris;

end

function [ris_coords, elements_per_ris] = splitRIS(original_center, total_elements, n)
    % 分割一个大RIS为n个小RIS
    % original_center: 大RIS的中心坐标 [x, y]
    % total_elements: 大RIS的总元素数量
    % n: 要分割的小RIS的数量
    % 返回小RIS的中心坐标数组和每个小RIS的元素数量
    
    % 每个小RIS的元素数量
    elements_per_ris = floor(total_elements / n);

    % 假设原来大RIS在x轴上的长度为1m
    length_x = 1;

    % 初始化小RIS的中心坐标数组
    ris_coords = zeros(n, 2);

    % 计算每个小RIS的中心坐标
    for i = 1:n
        ris_coords(i, :) = original_center + [(i - (n + 1) / 2) * length_x / n, 0];
    end
end

function [sum_outage_LOS, sum_outage_combined] = SA_Pout(Rx_coords, all_ris_coords, all_elements_per_ris, new_assignment)

% Simulation parameters
Lx = 15;
Ly = 10;
x_Tx = 8;
y_Tx = 3;
GHz = 18;
trials = 3000;

% RIS coordinates
% ris_coords = [0, 0; 0, 15; 20, 0; 20, 15];

s_light = 299792458;  % Speed of light
r = 0.1;  % Obstacle radius
lmbd = 0.5;  % Obstacle density
noise_variance = 1.658e-11;  % Noise variance
transmit_power = 0.1;  % Transmit power

% Antenna gains
Gt_dB = 20;  % Transmit antenna gain (dB)
Gr_dB = 10;  % Receive antenna gain (dB)
Gt = 10^(Gt_dB / 10);
Gr = 10^(Gr_dB / 10);
antenna_gain = Gr * Gt;
B = noise_variance / (transmit_power * antenna_gain);

% Frequency and attenuation settings
if GHz == 18
    z = 130;
    AO = 0.00006;
    f = 18e9;
elseif GHz == 26
    z = 200;
    AO = 0.00013;
    f = 26e9;
elseif GHz == 60
    z = 390;
    AO = 0.015;
    f = 60e9;
elseif GHz == 73
    z = 420;
    AO = 0.0075;
    f = 73e9;
end

wavelength = s_light / f;  % Signal wavelength
dx = wavelength / 2;  % RIS element size

num_Rx = size(Rx_coords, 1);

% Initialize outage probabilities
outage_LOS = zeros(num_Rx, 1);
outage_combined = zeros(num_Rx, 1);

for rx_idx = 1:num_Rx
    x_Rx = Rx_coords(rx_idx, 1);
    y_Rx = Rx_coords(rx_idx, 2);

    N = all_elements_per_ris(rx_idx, :);
    ris_idx = new_assignment(1, rx_idx);
    ris_coords = all_ris_coords(ris_idx, :);

    % Calculate LOS path loss
    [A_TxRx, B_TxRx, C_TxRx, ~] = line_equation_coefficients(x_Tx, x_Rx, y_Tx, y_Rx, 'single');
    distance_LOS = sqrt((x_Rx - x_Tx)^2 + (y_Rx - y_Tx)^2);
    path_loss_LOS = 20 * log10(4 * pi * f * distance_LOS / s_light) + AO * distance_LOS;

    % Initialize obstacle lengths
    obstacle_lengths_main_path = zeros(1, trials);
    obstacle_lengths_RIS_paths = zeros(N, trials, size(ris_coords, 1));

    % Simulate obstacle positions
    num_obstacles = poissrnd(lmbd * (Lx - 2 * r) * (Ly - 2 * r), [1, trials]);
    for trial = 1:trials
        x_obs = (Lx - 2 * r) * rand(1, num_obstacles(trial)) + r;
        y_obs = (Ly - 2 * r) * rand(1, num_obstacles(trial)) + r;

        % Calculate distances to obstructions
        distances_to_obstructions = dist_to_obstructions(x_obs, y_obs, A_TxRx, B_TxRx, C_TxRx);
        blocking_rays = blocked_rays_main_paths(x_obs, y_obs, x_Tx, x_Rx, y_Tx, y_Rx, r, distances_to_obstructions);

        obstacle_positions = r^2 - distances_to_obstructions.^2;
        obstacle_positions = obstacle_positions .* blocking_rays;
        obstacle_lengths_main_path(trial) = sum(2 * sqrt(obstacle_positions), 2);


        for idx = 1:size(ris_coords, 1)
            ris_center_x = ris_coords(idx, 1);
            ris_center_y = ris_coords(idx, 2);

            % Define RIS element positions
            x_RIS = dx * (linspace(1 - N / 2, N / 2, N) - 0.50) + ris_center_x;
            y_RIS = zeros(1, N);
            if ris_center_y ~= 0
                y_RIS = ones(1, N) * ris_center_y;
            end

            % Calculate paths
            [A_TxRIS, B_TxRIS, C_TxRIS, ~] = line_equation_coefficients(x_Tx, x_RIS, y_Tx, y_RIS, 'multiple');
            [A_RISRx, B_RISRx, C_RISRx, ~] = line_equation_coefficients(x_RIS, x_Rx, y_RIS, y_Rx, 'multiple');

            x_obs_repeated = repmat(x_obs', 1, N);
            y_obs_repeated = repmat(y_obs', 1, N);

            d2o_TxRIS = abs(repmat(A_TxRIS, num_obstacles(trial), 1) .* x_obs_repeated + repmat(B_TxRIS, num_obstacles(trial), 1) .* y_obs_repeated + repmat(C_TxRIS, num_obstacles(trial), 1)) ./ repmat(sqrt(A_TxRIS.^2 + B_TxRIS.^2), num_obstacles(trial), 1);
            d2o_RISRx = abs(repmat(A_RISRx, num_obstacles(trial), 1) .* x_obs_repeated + repmat(B_RISRx, num_obstacles(trial), 1) .* y_obs_repeated + repmat(C_RISRx, num_obstacles(trial), 1)) ./ repmat(sqrt(A_RISRx.^2 + B_RISRx.^2), num_obstacles(trial), 1);

            % Check for blocked rays between Tx-RIS and RIS-Rx
            blocking_rays_TxRIS = ((x_obs_repeated >= min(repmat(x_Tx, num_obstacles(trial), N), x_RIS) - r) & (x_obs_repeated <= max(repmat(x_Tx, num_obstacles(trial), N), x_RIS) + r) & ...
                                   (y_obs_repeated >= min(repmat(y_Tx, num_obstacles(trial), N), y_RIS) - r) & (y_obs_repeated <= max(repmat(y_Tx, num_obstacles(trial), N), y_RIS) + r) & ...
                                   (d2o_TxRIS < r));
            blocking_rays_RISRx = ((x_obs_repeated >= min(repmat(x_Rx, num_obstacles(trial), N), x_RIS) - r) & (x_obs_repeated <= max(repmat(x_Rx, num_obstacles(trial), N), x_RIS) + r) & ...
                                   (y_obs_repeated >= min(repmat(y_Rx, num_obstacles(trial), N), y_RIS) - r) & (y_obs_repeated <= max(repmat(y_Rx, num_obstacles(trial), N), y_RIS) + r) & ...
                                   (d2o_RISRx < r));

            % Calculate obstacle lengths
            obstacle_positions_TxRIS = 2 * sqrt(r^2 - d2o_TxRIS.^2) .* blocking_rays_TxRIS;
            obstacle_positions_RISRx = 2 * sqrt(r^2 - d2o_RISRx.^2) .* blocking_rays_RISRx;

            obstacle_lengths_RIS_paths(:, trial, idx) = sum(obstacle_positions_TxRIS + obstacle_positions_RISRx, 1);
        end
    end

    % Calculate path loss with obstacles
    path_loss_obstacles = 10.^(0.1 * (path_loss_LOS + z * obstacle_lengths_main_path));

    % Calculate path loss for RIS paths
    path_loss_RIS_list = [];
    for idx = 1:size(ris_coords, 1)
        ris_center_x = ris_coords(idx, 1);
        ris_center_y = ris_coords(idx, 2);

        x_RIS = dx * (linspace(1 - N / 2, N / 2, N) - 0.50) + ris_center_x;
        y_RIS = zeros(1, N);
        if ris_center_y ~= 0
            y_RIS = ones(1, N) * ris_center_y;
        end

        [~, ~, ~, d_TxRIS] = line_equation_coefficients(x_Tx, x_RIS, y_Tx, y_RIS, 'multiple');
        [~, ~, ~, d_RISRx] = line_equation_coefficients(x_RIS, x_Rx, y_RIS, y_Rx, 'multiple');

        Fcombine = fcombine('bottom', 'Conventional array', x_RIS, y_RIS, x_Tx, y_Tx, x_Rx, y_Rx, 1, (Gt/2)-1, (Gr/2)-1);
        Dair = sqrt(10.^(AO * (d_TxRIS + d_RISRx - 2 * r) / 10));
        num1 = repmat(sqrt(reshape(Fcombine,[],1)), 1, trials);
        LkRIS = 20 * log10(4 * pi / (dx^2));

        D_obs = sqrt(10.^(z * obstacle_lengths_RIS_paths(:, :, idx) / 10));
        den1 = (reshape(Dair, [], 1) .* D_obs .* reshape(d_TxRIS, [], 1) .* reshape(d_RISRx, [], 1));
        RIS_term = 20 * log10(abs(sum(num1 ./ den1, 1)));

        path_loss_RIS = zeros(2, trials);
        path_loss_RIS(1, :) = 10.^(0.1 * (path_loss_LOS + z * obstacle_lengths_main_path));
        path_loss_RIS(2, :) = 10.^(0.1 * (LkRIS - RIS_term));
        path_loss_RIS_list = min(path_loss_RIS);
    end

    path_loss_RIS_combined = path_loss_RIS_list;
    combined_path_loss = min(path_loss_obstacles, path_loss_RIS_combined);

    % Calculate outage probability
    sqr_h = exprnd(1, 1, trials);
    gamma_th = -10;
    len_gamma = length(gamma_th);

    compare_LOS = (repmat(sqr_h, len_gamma, 1) ./ repmat(path_loss_obstacles, len_gamma, 1))' <= B * 10.^(0.1 * gamma_th);
    outage_LOS(rx_idx) = sum(compare_LOS, 1) * (1 / trials);
    % avg_outage_LOS = mean2(outage_LOS);
    sum_outage_LOS = sum(outage_LOS);
    % % Confidence interval of LoS path
    % alpha_th = 0.05;
    % los_succ_count = sum(compare_LOS); % the number of success
    % los_z = norminv(1-alpha_th/2);  % the critical value of the normal distribution
    % p_hat_los = los_succ_count / trials; % the ratio of success
    % los_se = sqrt(p_hat_los * (1-p_hat_los) / trials); % Caculation the standrad error
    % ci_lower_los = p_hat_los - los_z * los_se;
    % ci_upper_los = p_hat_los + los_z * los_se;
    % CI_LOS = [ci_lower_los, ci_upper_los];

    compare_combined = (repmat(sqr_h, len_gamma, 1) ./ repmat(combined_path_loss, len_gamma, 1))' <= B * 10.^(0.1 * gamma_th);
    outage_combined(rx_idx) = sum(compare_combined, 1) * (1 / trials);
    % avg_outage_combined = mean2(outage_combined);
    sum_outage_combined = sum(outage_combined);
    % % Confidence interval of Combine path
    % comb_succ_count = sum(compare_combined); % the number of success
    % p_hat_comb = comb_succ_count / trials; % the ratio of success
    % comb_se = sqrt(p_hat_comb * (1-p_hat_comb) / trials); % Caculation the standrad error
    % ci_lower_comb = p_hat_comb - los_z * comb_se;
    % ci_upper_comb = p_hat_comb + los_z * comb_se;
    % CI_combined = [ci_lower_comb, ci_upper_comb];

    % fprintf('Receiver %d - LOS Outage Probability: %.4f, Combined Outage Probability: %.4f\n', rx_idx, outage_LOS(rx_idx), outage_combined(rx_idx));

end

function [A, B, C, d] = line_equation_coefficients(x_start, x_end, y_start, y_end, type_d)
    if strcmp(type_d, 'single')
        if (x_end - x_start) ~= 0
            m = (y_end - y_start) / (x_end - x_start);
        else
            m = 1e12;
        end
    elseif strcmp(type_d, 'multiple')
        m = (y_end - y_start) ./ (x_end - x_start);
        m(isinf(m)) = 1e12;
    end
    c = y_end - m .* x_end;
    B = ones(size(m));
    A = -m;
    C = -c;
    d = sqrt((x_end - x_start).^2 + (y_end - y_start).^2);
end

function distance = dist_to_obstructions(x_1, y_1, A, B, C)
    num = abs(A * x_1 + B * y_1 + C);
    den = sqrt(A.^2 + B.^2);
    distance = num / den;
end

function D3 = blocked_rays_main_paths(x, y, x1, x2, y1, y2, r, DS)
    minx = min([x1, x2]);
    maxx = max([x1, x2]);
    miny = min([y1, y2]);
    maxy = max([y1, y2]);
    D3 = ((x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy) & (DS <= r));
end

function Fcombine = fcombine(wall, RIS_type, x_RIS, y_RIS, x_Tx, y_Tx, x_Rx, y_Rx, alphaRIS, alphaTx, alphaRx)
    if strcmp(wall, 'bottom') || strcmp(wall, 'upper')
        x_center = abs(x_RIS(end) + x_RIS(1) + 2 * (x_RIS(2) - x_RIS(1))) / 2;
        if strcmp(wall, 'bottom')
            y_center = 0;
        else
            y_center = y_RIS(1);
        end
    else
        if strcmp(wall, 'left')
            x_center = 0;
        else
            x_center = x_RIS(1);
        end
        y_center = abs(y_RIS(end) + y_RIS(1) + 2 * (y_RIS(2) - y_RIS(1))) / 2;
    end

    theta1 = atan(abs(x_Tx - x_center) ./ abs(y_Tx - y_center)) - atan(abs(x_Tx - x_RIS) ./ abs(y_Tx - y_RIS));
    theta2 = atan(abs(x_Tx - x_RIS) ./ abs(y_Tx - y_RIS));
    theta3 = atan(abs(x_Rx - x_RIS) ./ abs(y_Rx - y_RIS));
    theta4 = atan(abs(x_Rx - x_center) ./ abs(y_Rx - y_center)) - atan(abs(x_Rx - x_RIS) ./ abs(y_Rx - y_RIS));

    F1 = cos(theta1).^alphaTx;
    F2 = cos(theta2).^alphaRIS;
    F3 = cos(theta3).^alphaRIS;
    F4 = cos(theta4).^alphaRx;

    if strcmp(RIS_type, 'Intelligent array')
        Fcombine = F1 .* F4;
    else
        Fcombine = F1 .* F2 .* F3 .* F4;
    end
end

end
