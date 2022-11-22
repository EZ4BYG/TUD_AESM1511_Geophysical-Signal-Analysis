%% Information
% Title: AESM1511-Assignment4.1
% Time: 29-10-2022
% Name: Boyu Gao + Miaoyang Yuan + Yongting Zhu
% Student Number: 5794676 + 5732476 + 5846528
% E-mail: B.Gao@student.tudelft.nl + yuanmy2318@gmail.com + Y.Zhu-44@student.tudelft.nl

%% Assignment 4.1
% Title: Fourier series for Box function and Gibbs phenomenon at discontinuities
% 1) Plot the box function and determine if it satisfies the 3 Dirichlet conditions
% 2) Use a finite Fourier series to approximate the box function 
% 3) Compute the energy in the difference and describe the Gibbs phenomenon  

%% Global parameters
% 1) Physical parameters
% T1                   the half of the time interval 
% T0_1                 the fundamental period of the box function: T0_1 = 4*T1
% T0_2                 the fundamental period of the box function: T0_2 = 6*T1
% w0_1                 the fundamental angular frequency of the box function: w0_1 = 2*pi/T0_1
% w0_2                 the fundamental angular frequency of the box function: w0_2 = 2*pi/T0_2
% a0_1                 the coefficient of the 0-th harmonic component: a0_1 = 2*T1/T0_1  
% a0_2                 the coefficient of the 0-th harmonic component: a0_2 = 2*T1/T0_2 
% num_points           number of points in the time interval 
% N                    different orders of the Fourier series
% 2) Figure parameters
% fontsize             the font size of title + x-axis + y-axis
% linewidth            the linewidth of x-axis + y-axis


%% Clean sheet
clear variables  % remove items from workspace, freeing up system memory
close all        % close any open figures
clc              % and clear the command window


%% Global parameters
% 1) Physical paramters
T1 = 1;
T0_1 = 4*T1;
T0_2 = 6*T1;
w0_1 = 2*pi / T0_1;
w0_2 = 2*pi / T0_2;
a0_1 = 2*T1 / T0_1;
a0_2 = 2*T1 / T0_2;
num_points = 201;
N = [1, 3, 7, 21, 47, 79];
% 2) Figure parameters
fontsize = 14;
linewidth = 2; 


%% Task 1: Plot the box function when T0 = 4*T1
% Predefine parameters
task1_t = linspace(-T0_1/2, T0_1/2, num_points);     % time: 1 x 201
task1_B = zeros(1, num_points);                      % Box: 1 x 201
% The Box function 
task1_B( abs(task1_t) <= T1 ) = 1;
% Plot
Figure1 = figure();       % open a new figure 
plot(task1_t, task1_B, 'Linewidth', linewidth);   
xlabel('time [s]', 'FontSize', fontsize);
ylabel('Amplitude', 'FontSize', fontsize);
title('T1. The exact box function (T_0 = 4T_1)', 'FontSize', fontsize);

% Displays
disp('Task1 quesiton: Does this function satisfy the three Dirichlet conditions?');
disp('Answer: Yes.');
disp('The exact box function is obviously absolutely integrable over the fundamental period T0,');
disp('and in each period, this function has finite max, min, discontinuous points.'); 
disp('---');


%% Task 2: Fourier series approximation when T0 = 4*T1
% Predefine parameters
task2_Ba = zeros(length(N), num_points);   % all kinds (6) of approximation: 6 x 201
task2_Energy = zeros(1, length(N));        % all kinds (6) of the energy in the difference: 1 x 6
task2_t = task1_t;                         % time: 1 x 201
num = 0;                                   % Counter
Figure2 = figure();                        % open a new figure 
sgtitle('T2-1. Fourier series of different orders approximation (T_0 = 4T_1)');  % total title
% Dual Loop
for N_tmp = N
    num = num + 1;                         % counter + 1
    an_tmp = zeros(1, N_tmp);              % the coefficients at the current N
    % Note: FS_tmp  —— The Fourier series has N_tmp items
    FS_tmp = zeros(N_tmp, num_points);     % Each row is the value of each item at all time points: N_tmp x 201
    % 1) Calculate the Fourier series
    for n_tmp = 1:N_tmp
        an_tmp(1, n_tmp) = sin(n_tmp*w0_1*T1) / (n_tmp*pi);                 % each coefficient
        FS_tmp(n_tmp, :) = 2 * an_tmp(1, n_tmp) * cos(n_tmp*w0_1*task2_t);  % each item of the Fourier series (all time)
    end
    % 2) Store each approximated result: sum of each col
    % Note: task2_Ba(num, :) means the sum of Fourier series at each point in time
    [row, col] = size(FS_tmp);
    if row == 1
        task2_Ba(num, :) = FS_tmp + a0_1;      
    else
        task2_Ba(num, :) = sum(FS_tmp) + a0_1;
    end
    % 3) Plot each approximation: 6 in total
    subplot(2, 3, num)
    plot(task1_t, task1_B, 'Linewidth', linewidth);           % the exact box function
    hold on                                                   % plot in one subfigure
    plot(task2_t, task2_Ba(num,:), 'Linewidth', linewidth);   % the approximatied box function
    % Mark the maxs and mins
    hold on
    max_y = max(task2_Ba(num,:));
    max_x = task2_t( find( task2_Ba(num,:) == max(task2_Ba(num,:)) ) );
    min_y = min(task2_Ba(num,:));
    min_x = task2_t( find( task2_Ba(num,:) == min(task2_Ba(num,:)) ) );
    plot(max_x, ones(1,length(max_x))*max_y, '*g', 'Linewidth', linewidth);  % maybe more than 1 max point
    hold on 
    plot(min_x, ones(1,length(min_x))*min_y, 'xb', 'Linewidth', linewidth);  % maybe more than 1 min point
    % Legend + Label + Subtitle 
    legend({'Exact', 'Approximated', 'Maxs', 'Mins'}, 'Location','south'); 
    xlabel('time [s]', 'FontSize', fontsize);
    ylabel('Amplitude', 'FontSize', fontsize);
    title(['N = ' num2str(N_tmp) ...
           ', B(+/-T_1) = ' num2str( task2_Ba(1, find(task2_t == -1)) ) ...
           ', B_{max}(t) = ' num2str( max(task2_Ba(num, :)) ) ]);     % each subfigure's title
end
% Calculate and plot the energy in difference 
for n_tmp = 1:length(N)
    square_error_tmp = (task2_Ba(n_tmp,:) - task1_B).^2;
    task2_Energy(1, n_tmp) = sum( square_error_tmp );
end
% Plot
Figure3 = figure();
plot(N, task2_Energy, 'Linewidth', linewidth);
xlabel('N', 'FontSize', fontsize);
ylabel('The energy in the difference', 'FontSize', fontsize);
title('T2-2. The energy in the difference when N changes (T_0 = 4T_1)', 'FontSize', fontsize);

% Displays
disp('Task2 question1: Does the result at the locations t = ± T1 change with N?');
disp('Answer: No. The values at t = ± T1 are always 0.5.');
disp('Task2 question2: Does the peak value change with N?');
disp('And explain what happens to the energy in the difference when N changes.');
disp('Answer: Yes. The peak value always changes with N. The bigger N, the smaller peak value.');
disp('Given the box function satisfies the three Dirichlet conditions, the Fourier series of this function must converge,');
disp('so as N increases, the approximating accuracy improves, the energy in the difference decreases, ');
disp('and the Gibbs phenomenon weakens.');
disp('---');


%% Task 3: Fourier series approximation when T0 = 6*T1
% Difference from Task 2: T0_1 → T0_2, w0_1 → w0_2, a0_1 → a0_2 
% Predefine parameters
task3_t = linspace(-T0_2/2, T0_2/2, num_points);    % T0_1 → T0_2 √
task3_B = zeros(1, num_points);
task3_B( abs(task3_t) <= T1 ) = 1;
task3_Ba = zeros(length(N), num_points);            % all kinds (6) of approximation: 6 x 201
task3_Energy = zeros(1, length(N));                 % all kinds (6) of the energy in the difference: 1 x 6
num = 0;                                            % Counter
Figure4 = figure();                                 % open a new figure 
sgtitle('T3-1. Fourier series of different orders approximation (T_0 = 6T_1)');  % total title
% Dual Loop
for N_tmp = N
    num = num + 1;                         % counter + 1
    an_tmp = zeros(1, N_tmp);              % the coefficients at the current N
    % Note: FS_tmp  —— The Fourier series has N_tmp items
    FS_tmp = zeros(N_tmp, num_points);     % Each row is the value of each item at all time points: N_tmp x 201
    % 1) Calculate the Fourier series
    for n_tmp = 1:N_tmp
        an_tmp(1, n_tmp) = sin(n_tmp*w0_2*T1) / (n_tmp*pi);                 % w0_1 → w0_2 √
        FS_tmp(n_tmp, :) = 2 * an_tmp(1, n_tmp) * cos(n_tmp*w0_2*task3_t);  % w0_1 → w0_2 √
    end
    % 2) Store each approximated result: sum of each col
    % Note: task2_Ba(num, :) means the sum of Fourier series at each point in time
    [row, col] = size(FS_tmp);
    if row == 1
        task3_Ba(num, :) = FS_tmp + a0_2;        % a0_1 → a0_2   
    else
        task3_Ba(num, :) = sum(FS_tmp) + a0_2;   % a0_1 → a0_2
    end
    % 3) Plot each approximation: 6 in total
    subplot(2, 3, num)
    plot(task3_t, task3_B, 'Linewidth', linewidth);           % the exact box function
    hold on                                                   % plot in one subfigure
    plot(task3_t, task3_Ba(num,:), 'Linewidth', linewidth);   % the approximatied box function
    % Mark the maxs and mins
    hold on
    max_y = max(task3_Ba(num,:));
    max_x = task3_t( find( task3_Ba(num,:) == max(task3_Ba(num,:)) ) );
    min_y = min(task3_Ba(num,:));
    min_x = task3_t( find( task3_Ba(num,:) == min(task3_Ba(num,:)) ) );
    plot(max_x, ones(1,length(max_x))*max_y, '*g', 'Linewidth', linewidth);  % maybe more than 1 max point
    hold on 
    plot(min_x, ones(1,length(min_x))*min_y, 'xb', 'Linewidth', linewidth);  % maybe more than 1 min point
    % Legend + Label + Subtitle 
    legend({'Exact', 'Approximated', 'Maxs', 'Mins'}, 'Location','south'); 
    xlabel('time [s]', 'FontSize', fontsize);
    ylabel('Amplitude', 'FontSize', fontsize);
    title(['N = ' num2str(N_tmp) ...
           ', B_{min}(t) = ' num2str( min(task3_Ba(num, :)) ) ...
           ', B_{max}(t) = ' num2str( max(task3_Ba(num, :)) ) ]);     % each subfigure's title
end
% Calculate and plot the energy in difference 
for n_tmp = 1:length(N)
    square_error_tmp = (task3_Ba(n_tmp,:) - task3_B).^2;
    task3_Energy(1, n_tmp) = sum( square_error_tmp );
end
% Plot
Figure5 = figure();
plot(N, task3_Energy, 'Linewidth', linewidth);
xlabel('N', 'FontSize', fontsize);
ylabel('The energy in the difference', 'FontSize', fontsize);
title('T3-2. The energy in the difference when N changes (T_0 = 6T_1)', 'FontSize', fontsize);

% Displays
disp('Task3 question: What happens when we take T0 = 6T1?');
disp('Answer: First, when we take T0 = 6T1, the Fourier series still converges.'); 
disp('Second, the time interval ([-T1,T1]) is smaller relative to the new period (T0=6T1),');
disp('so more high-frequency items (larger N) are needed to improve the approximation accuracy.');
disp('Compared with the T2.2, T-3.2 has a larger energy difference at the beginning,');
disp('but the energy differences decrease rapidly as N increases.');
disp('---')
