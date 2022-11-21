%% Information
% Title: AESM1511-Assignment3.1
% Time: 22-10-2022
% Name: Boyu Gao + Miaoyang Yuan
% Student Number: 5794676 + 5732476
% E-mail: B.Gao@student.tudelft.nl + yuanmy2318@gmail.com

%% Assignment 3.1
% Title: Filter surface waves in Frequency-wavenumber (f-k) domain 
% 1) Write a function file to read the original data and visualize it 
%    (scale the colorbar to 0.01 of the max and min value).
% 2) Use fft() and ifft() to transform the common-source gather from 
%    time domain to frequency-wavenumber domain. And explain why. 
% 3) Filter out the surface-wavenumber energy in frequency-wavenumber domain
%    using f-k filtering —— estimate the velocity of surface waves (a cone in the f-k plane)
%    , visualize the data in f-k domain and mute the energy we don't want.
% 4) Tranform the filtered data back to time-space domain. 

%% Global parameters
% 1) File parameters:
% fpath            the filepath of the original binary seimic data file 
% dat_format       the precision of the the binary file is float32
% mac_format       the machine format of the the binary file is IEEE floating point with little-endian byte ordering
% 2) Data parameters:
% traces           total numver of receivers
% time_samples     total samples collected by each receiver
% time_sampling    time interval for collecting data (used in fft)
% sampling_freq    = 1 / time_sampling
% num_freq         number of frequencies 
% freq_sampling    = sampling_freq / num_freq (used in ifft)
% space_sampling   space interval for receivers   (used in fftshift)
% receiver_start   the start space position (510m) of receivers 
% receiver_end     the end space position (1510m) of receivers
% 3) Filter parameters:
% filter           low-cut filter: frequencies <= 200 Hz
% vs               the estimated velocity of surface waves
% vr               the estimated velocity of reflection waves
% dv               velocity difference between vr and vs
% 4) Figure parameters: 
% fontsize_large   the font size of title + x-axis + y-axis
% fontsize_small   the font size of colorbar axis
% Linewidth        the linewidth of x-axis + y-axis
% colorbar_scale   scale the colorbar to 0.01 of max / min values

%% Main and Function file lists:
% 1) This file is the main function. 
% 2) ReadData_Func.m    Read binary file: input file + data parameters, return 2D array


%% Clean sheet
clear variables  % remove items from workspace, freeing up system memory
close all        % close any open figures
clc              % and clear the command window


%% Global parameters
% 1) File parameters
fpath = 'refl_3layers_fp50_dx0p5_500_rvz_tapered.bin';
dat_format = 'float32';
mac_format = 'ieee-le'; 
% 2) Data parameters
traces = 401;   
time_samples = 1001;
time_sampling = 0.001;
sampling_freq = 1 / time_sampling;
num_freq = time_samples; 
freq_sampling = sampling_freq / num_freq;
space_sampling = 2.5;
receiver_start = 510;
receiver_end = 1510;
% 3) Filter parameters
filter = 200;  
vs = 460;                
vr = 1200;              
dv = vr - vs;            
% 4) Figure parameters
fontsize_large = 14;
fontsize_small = 10;
Linewidth = 2;           
colorbar_scale = 0.01;


%% Task 1: Obtain and visualize the Reflection data 
% Obtain data
data_refl_time = ReadData_Func(fpath, traces, time_samples, dat_format, mac_format);
% Visualize 
x1_range = [receiver_start, receiver_end];   % horizontal distance from 510m to 1510m
y1_range = [0, time_samples];                % vertical time samples 
% Plot
Figure1 = figure();       % open a new figure 
imagesc(x1_range, y1_range, data_refl_time); % visualize the time-domain 2D (1001 x 401) reflection data
colormap('gray');         % use grey colour scale 
clb1 = colorbar;          % show the colour bar (gray), use clb1 to control the colorbar
xlabel('Horizontal Distances (m)', 'FontSize', fontsize_large);
ylabel('Time Samples', 'FontSize', fontsize_large);
title('A3.1-T1 Reflection data in the time domain', 'FontSize', fontsize_large);
ylabel(clb1, 'Amplitude', 'FontSize', fontsize_small);          % add a title to colorbar (at y location)
set(gca, 'linewidth', Linewidth);      % set the axis lines with linewidth of 2
set(gca, 'XAxisLocation', 'top');      % set the x-axis location to top (default is below)
set(gca, 'XTick', [receiver_start:250:receiver_end]);    % set the x-axis begin with 510 and end with 1510! step is 250
% Note: caxis is a new func name in R2022a —— set the limits of the current colorbar axis  
caxis( [colorbar_scale * min(min(data_refl_time)), colorbar_scale * max(max(data_refl_time))] );   

% Displays
disp('The size of the original time-domain reflection data: 1001 rows, 401 columns.');
disp('Each row: data collected by all 401 receivers at any time.')
disp('Each col: 1001 samples collected by each receiver.');
disp('---')


%% Task 2: Tranform the t-domain data to f-k domain 
% Transform
data_refl_fk = fftshift( ifft(fft(data_refl_time, [], 1)*time_sampling, [], 2)*(space_sampling),2 );  % complex
task2_fk = data_refl_fk( 1:filter, : );      % only need frequencies < 200 Hz
% Visualize
x2_range = [-0.5*time_sampling*traces, 0.5*time_sampling*traces];   % x-axis: wavenumber
y2_range = [0, filter];                                             % y-axis: frequency
% Plot
Figure2 = figure();       % open a new figure 
imagesc(x2_range, y2_range, abs(task2_fk));  % visualize the f-k domain 2D (1001 x 401) complex data 
colormap('default');               % use normal colour scale 
clb2 = colorbar;                   % show the colour bar (normal), use clb2 to control the colorbar
xlabel('Wavenumbers (1/m)', 'FontSize', fontsize_large);
ylabel('Frequency (Hz)', 'FontSize', fontsize_large);
title('A3.1-T2 Reflection data in the f-k domain', 'FontSize', fontsize_large);
ylabel(clb2, 'Amplitude', 'FontSize', fontsize_small);              % add a title to colorbar (at y location)
set(gca, 'linewidth', Linewidth);  % set the axis lines with linewidth of 2
set(gca, 'XAxisLocation', 'top');  % set the x-axis location to top (default is below)

% Displays
disp('Task2 question: Explain how to perform the f-k transformation using fft() + ifft().');
disp('Answer: fft2 can be replaced by fft twice (once fft, once ifft) —— ');
disp('First, fft in the time direction will convert the data from time-space domain to frequency-space domain;');
disp('Second, fft (ifft) in the space direction will convert the data from frequency-space domain to frequency-wavenumber domain.');
disp('Note: To better visualize the data in f-k domain, whether fft2 or fft twice, we need fftshift.');
disp('---');


%% Task 3: Mute the surface energy in f-k domain
% Note: Use cubic bandpass filter 
task3_fk = task2_fk;
task3_fk_size = size(task3_fk);
k_xaxis = linspace(-0.5*time_sampling*traces, 0.5*time_sampling*traces, task3_fk_size(2)); % x-axis: wavenumber —— the 201th is 0 (middle point)
f_yaxis = linspace(0, filter, task3_fk_size(1));                                           % y-axis: frequency
% Filter: set the energy (accordingly velocity from 1200 to 425) linearly to zero
% First: set the velocity < vs to zero —— which is lower than the surface waves
for col_tmp = 1:task3_fk_size(2)      % Column: 401
    for row_tmp = 1:task3_fk_size(1)  % Row: 200
        v_tmp = abs( f_yaxis(row_tmp) / k_xaxis(col_tmp) );
        if v_tmp < vs
            task3_fk(row_tmp, col_tmp) = 0;
        end
    end
end
% Second: linearly mute from vr to vs
for col_tmp = 1:task3_fk_size(2)      % Column: 401
    for row_tmp = 1:task3_fk_size(1)  % Row: 200
        v_tmp = abs( f_yaxis(row_tmp) / k_xaxis(col_tmp) );
        if (v_tmp >= vs) && (v_tmp <= vr)
            task3_fk(row_tmp, col_tmp) = task3_fk(row_tmp, col_tmp) * (v_tmp - vs) / dv;   % cubic bandpass filter
        end
    end
end
% Visualize
x3_range = x2_range;   % x-axis: wavenumber
y3_range = y2_range;   % y-axis: frequency
% Plot
Figure3 = figure();    % open a new figure 
imagesc(x3_range, y3_range, abs(task3_fk));  % visualize the f-k domain 2D (1001 x 401) complex data 
colormap('default');               % use normal colour scale 
clb2 = colorbar;                   % show the colour bar (normal), use clb2 to control the colorbar
xlabel('Wavenumbers (1/m)', 'FontSize', fontsize_large);
ylabel('Frequency (Hz)', 'FontSize', fontsize_large);
title('A3.1-T3 Reflection data in the f-k domain (filtered)', 'FontSize', fontsize_large);
ylabel(clb2, 'Amplitude', 'FontSize', fontsize_small);              % add a title to colorbar (at y location)
set(gca, 'linewidth', Linewidth);  % set the axis lines with linewidth of 2
set(gca, 'XAxisLocation', 'top');  % set the x-axis location to top (default is below)

% Displays
disp('Some important informain in f-k domain before and after filtering:');
disp("1) The magnitude of surface waves' amplitude is 10^(-6), the magnitude of reflection waves' amplitude is 10^(-8).");
disp('2) Two distinct energy lines with the same starting point (0,0), ');
disp('   the higher energy line is surface waves, the lower one is the reflection waves —— ');
disp('   which is obvious, because the amplitude of surface waves is much larger (100 times) than reflection waves.');
disp('3) Besides effective waves (surface and reflection waves) on the right, there are Aliasings on the left.');
disp('4) The slope of the line is velocity: surface waves have lower velocities, reflection waves have higher velocities.');
disp('   In contrast: in time domain, slope = y/x = 1/v; in f-k domain, slope = y/x = v.');
disp("5) Filtering is according to different slopes: cubic (slope's absolute value) bandpass filter —— ");
disp("   mute the energy where its slope's absolute value ( abs(y/x) ) is lower than the minimum reflected wave velocity.");
disp('6) Cubic bandpass filtering cannot filter out Aliasings.');
disp('7) After filtering and transform back to time domain, Aliasings are retained and appear in the upper right corner.');
disp('---');


%% Task4: Transform back to time domain
task4_fk = task3_fk;
% Transform 
task4_t = fft(ifft(task4_fk, [], 1), [], 2); % complex. I don't know why the value is complex.
% Visualize 
new_colorbar_scale = 0.0001;  % local paremter: 0.0001 better than 0.01 for visualizing the filtered time-domain data
x4_range = x1_range;          % horizontal distance from 510m to 1510m
y4_range = y1_range;          % vertical time samples 
% Plot
Figure4 = figure();      % open a new figure 
imagesc(x4_range, y4_range, real(task4_t));  % visualize the time-domain 2D (1001 x 401) reflection data (filtered)
colormap('gray');        % use grey colour scale 
clb4 = colorbar;         % show the colour bar (gray), use clb1 to control the colorbar
xlabel('Horizontal Distances (m)', 'FontSize', fontsize_large);
ylabel('Time Samples', 'FontSize', fontsize_large);
title('A3.1-T4 Reflection data in the time domain (filtered)', 'FontSize', fontsize_large);
ylabel(clb4, 'Amplitude', 'FontSize', fontsize_small);          % add a title to colorbar (at y location)
set(gca, 'linewidth', Linewidth);      % set the axis lines with linewidth of 2
set(gca, 'XAxisLocation', 'top');      % set the x-axis location to top (default is below)
set(gca, 'XTick', [receiver_start:250:receiver_end]);    % set the x-axis begin with 510 and end with 1510! step is 250
% Note: caxis is a new func name in R2022a —— set the limits of the current colorbar axis  
caxis( [new_colorbar_scale * min(min(data_refl_time)), new_colorbar_scale * max(max(data_refl_time))] ); 

% Displays
disp('Task4 question: What do you observe about the quality of the surface-wave events after filtering?');
disp('Answer: After testing, vs = 460 (m/s) is the best value for filtering out surface waves in this case.');
disp('1) In f-k domain: the energy line of surface waves is thrown out.');
disp('2) In time domain: the surfaces waves are almost completely deleted.');
disp('In sum, we can filter out surface waves perfectly in the f-k domain using cubic bandpass filter.');
disp('---');



