clear
clc

% Define positions
Tx = [8, 3];
RIS = [12, 0; 0, 6; 15, 5];
load("outage.mat");
load("user_trajectories.mat");


% Number of users and time slots
numUsers = 6;
S = [1, 10, 20, 30];
numTimeSlots = size(S, 2);

% Colors for each user
userColors = lines(numUsers);  % Generates a colormap with distinct colors

% Create figure with subplots
figure;
for t = 1:numTimeSlots
    subplot(2, 2, t);
    hold on;
    title(['Time Slot ', num2str(S(t))]);
    xlim([0, 15]);
    ylim([0, 10]);
    % xlabel('X Position');
    % ylabel('Y Position');
    grid on;
    
    % Plot Tx position
    plot(Tx(1), Tx(2), '^','color', 'black', 'MarkerSize', 7, 'MarkerFaceColor','black', 'LineWidth', 1.3, 'DisplayName', 'Tx');
    
    % Plot RIS positions
    plot(RIS(:, 1), RIS(:, 2), 'square', 'color', [18 133 66]/255, 'MarkerFaceColor',[18 133 66]/255, 'MarkerSize', 7, 'LineWidth', 1.3, 'DisplayName', 'RIS');
    
    % Plot users and connections to RIS
    for u = 1:numUsers
        userPos = trajectories{u}(S(t), 1:2);
        risAssoc = results{t}.assignment(1, u);
        
        % Plot user position with different colors
        if t == 1 && u==1  % Only add legend for the first time slot
            plot(userPos(1), userPos(2), 'o', 'Color', userColors(u, :), 'MarkerSize', 5, 'LineWidth', 1.3, 'DisplayName', 'Rx','HandleVisibility','on');
        else
            plot(userPos(1), userPos(2), 'o', 'Color', userColors(u, :), 'MarkerSize', 5, 'LineWidth', 1.3, 'HandleVisibility','off');
        end

        % Connect user to Tx
        if t == 1 && u==1
            plot([userPos(1), Tx(1)], [userPos(2), Tx(2)], 'color', '#9EC0E3', 'LineWidth', 1.3, 'DisplayName', 'LoS links', 'HandleVisibility','on');
        else
            plot([userPos(1), Tx(1)], [userPos(2), Tx(2)], 'color', '#9EC0E3', 'LineWidth', 1.3, 'HandleVisibility','off');
        end
        
        
        % Connect user to RIS if associated
        if risAssoc > 0 && t==1 && u==1
            plot([userPos(1), RIS(risAssoc, 1)], [userPos(2), RIS(risAssoc, 2)], '--','color', '#F17F7E', 'LineWidth', 1.3, 'DisplayName', 'Virtual links','HandleVisibility','on');
        else
            plot([userPos(1), RIS(risAssoc, 1)], [userPos(2), RIS(risAssoc, 2)], '--', 'color', '#F17F7E', 'LineWidth', 1.3, 'HandleVisibility','off');
        end
    end
    
    if t == 1
        h=legend('show', 'Location', 'BestOutside');
        set(h, 'Orientation', 'horizon', 'Box', 'on');

    end
end
