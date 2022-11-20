%% Information
% Title: AESM1511-Assignment2
% Time: 08-10-2022
% Name: Boyu Gao + Miaoyang Yuan
% Student Number: 5794676 + 5732476
% E-mail: B.Gao@student.tudelft.nl + yuanmy2318@gmail.com
%
%% Assignment 2
% Title: 2D Fourier Transformation and Filtering
% 1) Program to read (fopen/fread) and visualize (imagesc) the reflection seismic data;
% 2) use Fast Fourier Transform (fft) to transform the time-domain data to
% the frequency-space domain;
% 3) filter low frequencies in the frequency-space domain, and then use
% Inverse Fast Fourier Transform (ifft) to transform the filtered data back
% to the time-space domain;
% 4) transform the data from the frequency-space domain to the
% frequency-wavenumber domain (fftshift).

% Task1 + Task4 = gray      Task2 + 5 + 6 = default

%% Global parameters 
% fpath            the filepath of the original binary seimic data file 
% traces           total receivers from 510m to 1510m spaced at 2.5m 
% time_samples     total samples collected by each receiver
% time_sampling    time interval for collecting data (used in fft)
% sampling_freq    = 1 / time_sampling
% num_freq         number of frequencies 
% freq_sampling    = sampling_freq / num_freq (used in ifft)
% space_sampling   space interval for receivers   (used in fftshift)
% filter           low-cut filter: frequencies <= 200 Hz


%% Clean sheet
clear variables  % remove items from workspace, freeing up system memory
close all        % close any open figures
clc              % and clear the command window


%% Global parameters
fpath = 'refl_3layers_fp50_dx0p5_500_rvz.bin'; 
traces = 401;   
time_samples = 1001;
time_sampling = 0.001;
sampling_freq = 1 / time_sampling;
num_freq = time_samples; 
freq_sampling = sampling_freq / num_freq;
space_sampling = 2.5;
filter = 200;      


%% Reflection data 
fid = fopen(fpath, 'r');  % open the file
temp = fread(fid, traces * time_samples, 'float32', 0, 'ieee-le');  % fread read the binary file
fclose(fid);              % close the file  
data_refl_time(:,:) = reshape(temp, time_samples, traces);  % get the time-domain data and reshape its size to 1001 x 401
% Displays
disp('The size of the time-domain reflection data: 1001 rows, 401 columns.');
disp('Each row: data collected by all 401 receivers at any time.')
disp('Each col: 1001 samples collected by each receiver.');
disp('---')


%% Task 1
% Local parameters
x1_range = [510, 1510];   % horizontal distance from 510m to 1510m
y1_range = [0, 1001];     % vertical time samples 

% Plot
Figure1 = figure();       % open a new figure 
imagesc(x1_range, y1_range, data_refl_time);  % visualize the time-domain 2D (1001 x 401) reflection data
colormap('gray');         % use grey colour scale 
clb1 = colorbar;          % show the colour bar (gray), use clb1 to control the colorbar
xlabel('Horizontal Distance (m)', 'FontSize', 14);
ylabel('Time Samples', 'FontSize', 14);
title('T1. Reflection data in the time-domain', 'FontSize', 14);
ylabel(clb1, 'Amplitude', 'FontSize', 10);          % add a title to colorbar (at y location)
set(gca, 'linewidth', 2);             % set the axis lines with linewidth of 2
set(gca, 'XAxisLocation', 'top');     % set the x-axis location to top (default is below)
set(gca, 'XTick', [510:250:1510]);    % set the x-axis begin with 510 and end with 1510!


%% Task 2 
% Local parameters
x2_range = [510, 1510];   % horizontal (x) range
y2_range = [0, filter];   % vertical (y) range
dim_fft = 1;              % an important parameter in fft 
% FFT: fft(data, [], dimension)
% dimension = 2 means that operates along the rows and returns the Fourier transform of each row.
% dimension = 1 means that operates along the columns and returns the Fourier transform of each column.
data_refl_freq = fft(data_refl_time, [], dim_fft)*(time_sampling);  % the return value is complex number.
data_refl_freq( filter+1:end, : ) = 0;        % Filter: cut off the positive frequencies > 200 Hz
task2_freq = data_refl_freq( 1:filter, : );   % only show frequencies < 200 Hz

% Plot
Figrue2 = figure();                               % open a new figure
imagesc( x2_range, y2_range, abs(task2_freq) );   % absolute value is the amplitude!!
colormap('default');               % use normal colour scale 
clb2 = colorbar;                   % show the colour bar (normal), use clb2 to control the colorbar
xlabel('Horizontal Distance (m)', 'FontSize', 14);
ylabel('Frequency (Hz)', 'FontSize', 14);
title('T2. Reflection data in the frequency-domain (≤ 200 Hz)', 'FontSize', 14);
ylabel(clb2, 'Amplitude', 'FontSize', 10);            % add a title to colorbar (at y location)
set(gca, 'linewidth', 2);          % set the axis lines with linewidth of 2
set(gca,'XAxisLocation','top');    % set the x-axis location to top (default is below)
set(gca,'XTick', [510:250:1510]);

% Displays
disp('Answer to the question in Task 2')
disp('Q: What happens when attempting to plot and why that happens?');
disp('A: An error occurs if we directly visualize the fft result, cuz the result is complex-valued,');
disp('we should use abs(fft_result) as the input of imagesc.');
disp("If we correctly visualize the fft result, we will find 'frequency symmetry'.");
disp('The reason is the fourier series are positive-negative symmetry.')
disp('---');


%% Task 3 
% Display
disp('Answer to the question in Task 2')
disp('Q: What events in the time-space domain does that energy correspond?');
disp('A: The image shows that the energy is concentrated in the first 175 receivers,');
disp('The direct surface wave has not yet travelled to further receivers.');
disp('---')


%% Task 4
% Local parameters
x4_range = [510, 1510];                % horizontal (x) range
y4_range = [0, filter];                % vertical (y) range
filter_low = 30;                       % (2)(3) set the frequencies < 30 Hz to zero
filter_high = 60;                      % don't change the frequencies > 60 Hz
task4_freq = data_refl_freq;           % the data (complex) needed to be filtered is the result of Task 2
dfilter = filter_high - filter_low;    % distance between the filter_high and filter_low
x_lim = 175;                           % we only filter the first 175 receivers (adjustable)
dim_ifft = 1;                          % an important parameter in ifft

% Filter: 1001 x 401
% Filter (1): set the frequencies < 60 Hz to zero
task4_freq_1 = task4_freq;
task4_freq_1(1:filter_high-1, 1:x_lim) = 0; 
% Filter (2): set the frequencies from 60 Hz to 30 Hz go linearly to zero
% First: set the frequencies < 30 Hz to zero
% Second: times by a factor that linearly approaches to 1 
task4_freq_2 = task4_freq;
task4_freq_2(1:filter_low-1, 1:x_lim) = 0;   % First
row_tmp = 0;
for freq = filter_low:filter_high
    row_tmp = row_tmp + 1;
    factor = row_tmp / (dfilter+1);
    task4_freq_2(freq, 1:x_lim) = task4_freq_2(freq, 1:x_lim) * factor;  % Second
end
% Filter (3): set the frequencies from 60 Hz to 30 Hz go sinusoidally to zero
% First:  set the frequencies < 30 Hz to zero
% Second: times by a factor that sinusoidally approaches to 1 
task4_freq_3 = task4_freq;
task4_freq_3(1:filter_low-1, 1:x_lim) = 0;  % First
row_tmp = 0;
for freq = filter_low:filter_high
    row_tmp = row_tmp + 1;
    factor = sin( row_tmp / (dfilter+1) * 4 / pi );        % T/4
    task4_freq_3(freq, 1:x_lim) = task4_freq_3(freq, 1:x_lim) * factor;  % Second
end
% IFFT: 1001 x 401
task4_time_1 = 2*real( num_freq * ifft(task4_freq_1, [], dim_ifft) * freq_sampling );
task4_time_2 = 2*real( num_freq * ifft(task4_freq_2, [], dim_ifft) * freq_sampling );
task4_time_3 = 2*real( num_freq * ifft(task4_freq_3, [], dim_ifft) * freq_sampling );
% Put data in cell for loop plotting 
% Time-domain: 1001 x 401; Frequency-domain: 200 x 401
task4_time_cell = { task4_time_1, task4_time_2, task4_time_3 };
task4_freq_cell = { task4_freq_1(1:filter,:), task4_freq_2(1:filter,:), task4_freq_3(1:filter,:) };

% Plot1: Frequency-domain (after FFT)
title1 = {'T4. Frequency-domain: Zero filtering', ...
          'T4. Frequency-domain: Linear filtering', ...
          'T4. Frequency-domain: Sin filtering'};
for f_tmp = 1:3
    Figrue4 = figure();                     % open a new figure
    imagesc(x4_range, y4_range, abs(task4_freq_cell{1,f_tmp}) );    % 200 x 401
    colormap('default');                    % use normal colour scale 
    clb4 = colorbar;                        % show the colour bar (normal), use clb2 to control the colorbar
    xlabel('Horizontal Distance (m)', 'FontSize', 14);
    ylabel('Frequency (Hz)', 'FontSize', 14);
    title( title1{1,f_tmp}, 'FontSize', 14);
    ylabel(clb4, 'Amplitude', 'FontSize', 10);    % add a title to colorbar (at y location)
    set(gca, 'linewidth', 2);               % set the axis lines with linewidth of 2
    set(gca, 'XAxisLocation', 'top')        % set the x-axis location to top (default is below)
    set(gca, 'XTick', 510:250:1510);
end
% Plot2: Time-domain (after IFFT)
title2 = {'T4. Time-domain: After zero filtering', ...
          'T4. Time-domain: After linear filtering', ...
          'T4. Time-domain: After sin filtering'};
for f_tmp = 1:3
    Figrue4 = figure();                     % open a new figure
    imagesc(x1_range, y1_range, task4_time_cell{1,f_tmp} );    % 1001 x 401
    colormap('gray');                       % use normal colour scale 
    clb4 = colorbar;                        % show the colour bar (normal), use clb2 to control the colorbar
    xlabel('Horizontal Distance (m)', 'FontSize', 14);
    ylabel('Time samples', 'FontSize', 14);
    title( title2{1,f_tmp}, 'FontSize', 14);
    ylabel(clb4, 'Amplitude', 'FontSize', 10);     % add a title to colorbar (at y location)
    set(gca, 'linewidth', 2);               % set the axis lines with linewidth of 2
    set(gca, 'XAxisLocation','top');        % set the x-axis location to top (default is below)
    set(gca, 'XTick', [510:250:1510]);
end

% Display
disp('Answer to the question in Task 4');
disp('Q1: Should this be done for all horizontal distances?');
disp('A1: No, we just need to do the low-cut filters for the first 175 receivers.');
disp('Q2: What events are preserved and what is their preservation quality? Which filter gave the clearest result?');
disp('A2: These 3 filters reduces the low noise in the original data. The direct surface wave is preserved.');
disp('Linear and sinusoidal filters achieve almost the same good preservation quality,');
disp('but the filter setting the frequencies lower than 60 Hz to zero is the worst.');
disp('---');
 

%% Task 5 and 6
% Local parameters
task5_freq = data_refl_freq;    % 1001 x 401
dim_fft_2 = 2;                  % different from task 2 
x5_range = [-0.5*time_sampling*traces, 0.5*time_sampling*traces];
y5_range = [0, filter];
trace_num = {'1st', '2nd', '4th', '8th'};
% FFTSHIFT: transform data from the frequency-space domain to the frequency-wavenumber domain
% Data size: 1001 x 401
task5_freq_wavenum1 = fftshift( fft(task5_freq(:,1:end), [], dim_fft_2) * space_sampling, dim_fft_2);
task5_freq_wavenum2 = fftshift( fft(task5_freq(:,1:2:end), [], dim_fft_2) * space_sampling, dim_fft_2);
task5_freq_wavenum4 = fftshift( fft(task5_freq(:,1:4:end), [], dim_fft_2) * space_sampling, dim_fft_2);
task5_freq_wavenum8 = fftshift( fft(task5_freq(:,1:8:end), [], dim_fft_2) * space_sampling, dim_fft_2);
% Put data in cell for loop plotting: 200 x 401, 200 x 201, 200 x 101, 200 x 51
task5_freq_wavenum = {task5_freq_wavenum1(1:filter, :), ...
                      task5_freq_wavenum2(1:filter, :), ...
                      task5_freq_wavenum4(1:filter, :), ...
                      task5_freq_wavenum8(1:filter, :)};

% Plot
Figure5 = figure(); 
for f_tmp = 1:4
    subplot(2, 2, f_tmp);
    imagesc( x5_range, y5_range, abs(task5_freq_wavenum{1,f_tmp}) );
    colormap('default');               % use normal colour scale 
    clb5 = colorbar;                   % show the colour bar (normal), use clb2 to control the colorbar
    xlabel('Wavenumbers (1/m)', 'FontSize', 14);
    ylabel('Frequency (Hz)', 'FontSize', 14);
    title( "T5. Frequecny wavenumber domain: Every " + trace_num{1,f_tmp} + ' tracks', 'FontSize', 14);
    ylabel(clb5, 'Amplitude', 'FontSize', 10);          % add a title to colorbar (at y location)
    set(gca, 'linewidth', 2);          % set the axis lines with linewidth of 2
    set(gca, 'XAxisLocation', 'top');  % set the x-axis location to top (default is below)
end

% Display
disp('Answer to the question in Task 5 and 6');
disp('Q1: What can be observed about the energy?');
disp('A1: In frequency-wavenumber domain (after ft2 + fftshift), the center of the image represents low-frequency energy,');
disp(' the four corners of the image represent the high-frequency energy,');
disp('In task 5, the main energy is concentrated at the relatively high frequencies.');
disp('Q2: Comparing the 4 images, explain what happens to the main energy.');
disp('A2: As the traces decrease (the sampling becomes more sparse), more energy stripes appear in the image,');
disp('the main energy becomes more concentrated at the relatively high frequencies.');
disp('---');

