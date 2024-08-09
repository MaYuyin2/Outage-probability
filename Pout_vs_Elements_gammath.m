clc;
clear;

% Simulation parameters
Lx = 15;
Ly = 10;
x_Tx = 8;
y_Tx = 3;
% x_Rx = 14;
% y_Rx = 10;
x_Rx = 7;
y_Rx = 6;
GHz = 18;
trials = 10000;

% RIS coordinates
ris_coords = [12,0];
% ris_coords = [12,0; 15,5];

% Different N values to simulate
N_values = [64, 256];

% Pre-allocate space for storing outage probabilities and confidence intervals
outage_LOS_store = zeros(length(N_values), 201);
CI_LOS_store = zeros(length(N_values), 2, 201);
outage_combined_store = zeros(length(N_values), 201);
CI_combined_store = zeros(length(N_values), 2, 201);

for i = 1:length(N_values)
    N = N_values(i);
    [outage_LOS_store(i, :), CI_LOS_store(i, :, :), outage_combined_store(i, :), CI_combined_store(i, :, :)] = ...
        simulate(Lx, Ly, x_Tx, y_Tx, x_Rx, y_Rx, ris_coords, N, GHz, trials);
end

% Plot results
figure;
axes2=axes('position',[0.2,0.2,0.74,0.7]); %这个是figure里面图的位置和大小，分别为离下边，左边的距离，还有图的高和宽
set(gca, "LooseInset", [0,0,0,0]);%消除白边
hold on;
colors = {[254 129 125]/255, [129 184 223]/255, [72 96 170]/255, [241 127 126]/255};  % Define colors for different N values

% Processing the data (CI_LOS and CI_combined). The figure can not be 
% exhibited because there are 0 and negativer numbers
epsilon = 1e-10; % 
CI_LOS_store(CI_LOS_store <= 0) = epsilon;
CI_combined_store(CI_combined_store <= 0) = epsilon;

for i = 1:length(N_values)
    % Plot LOS results
    fill([linspace(-100, 100, 201), fliplr(linspace(-100, 100, 201))], ...
        [squeeze(CI_LOS_store(i, 1, :))', fliplr(squeeze(CI_LOS_store(i, 2, :))')], ...
        colors{i}, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(linspace(-100, 100, 201), outage_LOS_store(i, :), '--', 'LineWidth', 1.3, ...
        'Color', colors{i}, 'DisplayName', sprintf('LOS (N = %d)', N_values(i)));

    % Plot combined results
    fill([linspace(-100, 100, 201), fliplr(linspace(-100, 100, 201))], ...
        [squeeze(CI_combined_store(i, 1, :))', fliplr(squeeze(CI_combined_store(i, 2, :))')], ...
        colors{i}, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(linspace(-100, 100, 201), outage_combined_store(i, :), '-', 'LineWidth', 1.3, ...
        'Color', colors{i}, 'DisplayName', sprintf('LOS + RIS (N = %d)', N_values(i)));
end

set(gca, 'YScale', 'log');
xlim([-60, 100]);
ylim([0.001, 1]);
grid on;
set(gca, "FontSize",10, "Fontname", "Times new roman"); %轴刻度标签的字体大小和名称
ylabel('Outage Probability', "FontSize",12, "Fontname", "Times new roman");
xlabel('\rho_{th} [dB]', "FontSize",12, "Fontname", "Times new roman");
% title(sprintf('f = %dGHz, No trials = %d', GHz, trials));
h = legend('show', "Fontname", "Times new roman", 'location','northwest', "FontSize",10, 'edgecolor', [1,1,1]);
set(h, 'box', 'off'); %设置legend背景色透明
hold off;

function [outage_LOS, CI_LOS, outage_combined, CI_combined] = simulate(Lx, Ly, x_Tx, y_Tx, x_Rx, y_Rx, ris_coords, N, GHz, trials)
    c = 299792458;  % Speed of light
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

    wavelength = c / f;  % Signal wavelength
    dx = wavelength / 2;  % RIS element size

    % Calculate LOS path loss
    [A_TxRx, B_TxRx, C_TxRx, d] = line_equation_coefficients(x_Tx, x_Rx, y_Tx, y_Rx, 'single');
    distance_LOS = sqrt((x_Rx - x_Tx)^2 + (y_Rx - y_Tx)^2);
    path_loss_LOS = 20 * log10(4 * pi * f * distance_LOS / c) + AO * distance_LOS;

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
            [A_TxRIS, B_TxRIS, C_TxRIS, d_TxRIS] = line_equation_coefficients(x_Tx, x_RIS, y_Tx, y_RIS, 'multiple');
            [A_RISRx, B_RISRx, C_RISRx, d_RISRx] = line_equation_coefficients(x_RIS, x_Rx, y_RIS, y_Rx, 'multiple');

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

        [A_TxRIS, B_TxRIS, C_TxRIS, d_TxRIS] = line_equation_coefficients(x_Tx, x_RIS, y_Tx, y_RIS, 'multiple');
        [A_RISRx, B_RISRx, C_RISRx, d_RISRx] = line_equation_coefficients(x_RIS, x_Rx, y_RIS, y_Rx, 'multiple');

        Fcombine = fcombine('bottom', 'Conventional array', x_RIS, y_RIS, x_Tx, y_Tx, x_Rx, y_Rx, 1, (Gt/2)-1, (Gr/2)-1);
        Dair = sqrt(10.^(AO * (d_TxRIS + d_RISRx - 2 * r) / 10));
        num = repmat(sqrt(reshape(Fcombine,[],1)), 1, trials);
        LkRIS = 20 * log10(4 * pi / (dx^2));

        D_obs = sqrt(10.^(z * obstacle_lengths_RIS_paths(:, :, idx) / 10));
        den = (reshape(Dair, [], 1) .* D_obs .* reshape(d_TxRIS, [], 1) .* reshape(d_RISRx, [], 1));
        RIS_term = 20 * log10(abs(sum(num ./ den, 1)));

        path_loss_RIS = zeros(2, trials);
        path_loss_RIS(1, :) = 10.^(0.1 * (path_loss_LOS + z * obstacle_lengths_main_path));
        path_loss_RIS(2, :) = 10.^(0.1 * (LkRIS - RIS_term));
        path_loss_RIS_list = min(path_loss_RIS);
    end

    % path_loss_RIS_matrix= cell2mat(path_loss_RIS_list'); % coverted the cell array into list
    path_loss_RIS_combined = path_loss_RIS_list;
    combined_path_loss = min(path_loss_obstacles, path_loss_RIS_combined);

    % Calculate outage probability
    sqr_h = exprnd(1, 1, trials);
    gamma_th = linspace(-100, 100, 201);
    len_gamma = length(gamma_th);

    compare_LOS = (repmat(sqr_h, len_gamma, 1) ./ repmat(path_loss_obstacles, len_gamma, 1))' <= B * 10.^(0.1 * gamma_th);
    outage_LOS = sum(compare_LOS, 1) * (1 / trials);
    % Confidence interval of LoS path
    alpha_th = 0.05;
    los_z = norminv(1-alpha_th/2);  % the critical value of the normal distribution
    for i=1:len_gamma
        los_succ_count = sum(compare_LOS(:, i)); % the number of success
        p_hat_los = los_succ_count / trials; % the ratio of success
        los_se = sqrt(p_hat_los * (1-p_hat_los) / trials); % Caculation the standrad error
        ci_lower_los = p_hat_los - los_z * los_se;
        ci_upper_los = p_hat_los + los_z * los_se;
        CI_LOS(1,i) = ci_lower_los;
        CI_LOS(2,i) = ci_upper_los;
    end

    compare_combined = (repmat(sqr_h, len_gamma, 1) ./ repmat(combined_path_loss, len_gamma, 1))' <= B * 10.^(0.1 * gamma_th);
    outage_combined = sum(compare_combined, 1) * (1 / trials);
    % Confidence interval of Combine path
    for i=1:len_gamma
        comb_succ_count = sum(compare_combined(:,i)); % the number of success
        p_hat_comb = comb_succ_count / trials; % the ratio of success
        comb_se = sqrt(p_hat_comb * (1-p_hat_comb) / trials); % Caculation the standrad error
        ci_lower_comb = p_hat_comb - los_z * comb_se;
        ci_upper_comb = p_hat_comb + los_z * comb_se;
        CI_combined(1,i) = ci_lower_comb;
        CI_combined(2,i) = ci_upper_comb;
    end
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
