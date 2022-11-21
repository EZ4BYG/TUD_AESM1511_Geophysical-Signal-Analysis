%% Information
% Title: AESM1511-Assignment3.2
% Time: 24-10-2022
% Name: Boyu Gao + Miaoyang Yuan
% Student Number: 5794676 + 5732476
% E-mail: B.Gao@student.tudelft.nl + yuanmy2318@gmail.com

%% Assignment 3.2
% Title: Interpolation and Integration of Diffusive Fields
% 1) Create a function file to generate the integrand.
% 2) Use the integrand function file to calculate a series of frequencies 
%    with a logarithmic spacing from f = 10^(-2.5) Hz to f = 10^(2.5) Hz 
%    spanning 10^5 (5 decades), using 2 points per decade.
% 3) Calculate the interpolated results and use normalized values to test 
%    the accuray of the interpolated results. 
% 4) Bisect and obtain new frequency points through loop, and then use
%    these frequencies to compute the integral. 
% 5) Compute the field for all linear frequencies by interpolating the 
%    field computed at the limitd irregular frequencies obtained in Task 5.

%% Global parameters
% 1) Physical parameters:
% sigma     conductivity 
% mu        magnetic permeability
% rho       the horizontal distance 
% h         the veritcal distance
% f_start   the starting vale of the frequency
% f_end     the end value of the frequency

%% Main and function file lists:
% 1) This file is the man function. 
% 2) Integrand_Func.m     Generate the integrand: input physical paramters, return the integrand


%% Clean sheet
clear variables  % remove items from workspace, freeing up system memory
close all        % close any open figures
clc              % and clear the command window


%% Global parameters
% 1) Physical parameters 
sigma = 1;                       
mu = 4 * pi * 10^(-7);         
rho = 1000;                    
h = 500;
f_start = 10^(-2.5);
f_end = 10^(2.5);


%% Task 1 and 2: Create a function file and use it to calculate the integration 
% Local parameters
task1_f = [0.1, 1, 10];                     % frequencies 
task1_omega = zeros(1, length(task1_f));    % Initial: angular frequency
task1_gamma = zeros(1, length(task1_f));    % Initial: gamma —— one parameter in J0 
task1_F1 = zeros(1, length(task1_f));       % Initial: the integration results using J0 Bessel function
% Calculate the J0 integration result: complex!
for f_tmp = 1:length(task1_f)
    task1_omega(f_tmp) = 2 * pi * task1_f(f_tmp);                        % angular frequency
    task1_gamma(f_tmp) = sqrt( 1i * task1_omega(f_tmp) * sigma * mu );   % gamma
    % Obtain the integrand
    Integrand = @(kappa) Integrand_Func(kappa, task1_gamma(f_tmp), rho, h);  % function handle of the integrand 
    task1_F1(f_tmp) = 1i * task1_omega(f_tmp) * quadgk(Integrand, 0, inf);   % the integration result
end
% Calculate the I0 + K0 integration result: complex!
I0 = besseli( 0, 0.5*task1_gamma*(sqrt((rho^2)+(h^2))-h) );   % the first kind of modified zeroth Bessel function 
K0 = besselk( 0, 0.5*task1_gamma*(sqrt((rho^2)+(h^2))+h) );   % the second kind of modified zeroth Bessel function
task1_F2 = 1i * task1_omega .* I0 .* K0;
% Calculate the errors between F1 and F2
task1_Error = abs( task1_F1 - task1_F2  );   
% Display
fprintf('Task 1 and 2:\n');
for f_tmp = 1:length(task1_f)
    fprintf('When f = %.1f: ', task1_f(f_tmp));
    fprintf('F1 = %f + %fi, ', task1_F1(f_tmp), task1_F1(f_tmp)/1i);
    fprintf('F2 = %f + %fi, ', task1_F2(f_tmp), task1_F2(f_tmp)/1i);   
    fprintf('Error = %e\n', task1_Error(f_tmp));
end
disp('The error between the two formulas is very small, so we can reasonablely believe the two formulas are the same.');
disp('---');


%% Task 3: Calculate the integration results over a series of frequencies 
% Local parameters
task3_f_start = f_start;      % the start frequency is 10^(-2.5) Hz
task3_f_end = f_end;          % the end frequency is 10^(2.5) Hz
task3_space = 10^(0.5);       % space
task3_samples = log10(task3_f_end/task3_f_start) / log10(task3_space) + 1;   % total samples 
task3_f = task3_f_start * task3_space .^ (0:task3_samples-1);                % generate geometric series
% Calculate the integration over all frequencies 
task3_omega = zeros(1, task3_samples);    % Initial: angular frequency
task3_gamma = zeros(1, task3_samples);    % Initial: gamma —— one parameter in J0 
task3_F1 = zeros(1, task3_samples);       % Initial: the integration results using J0 Bessel function
for f_tmp = 1:task3_samples
    task3_omega(f_tmp) = 2 * pi * task3_f(f_tmp);                        % angular frequency
    task3_gamma(f_tmp) = sqrt( 1i * task3_omega(f_tmp) * sigma * mu );   % gamma
    % Obtain the integrand
    Integrand = @(kappa) Integrand_Func(kappa, task3_gamma(f_tmp), rho, h);  % function handle of the integrand 
    task3_F1(f_tmp) = 1i * task3_omega(f_tmp) * quadgk(Integrand, 0, inf);   % the integration result
end 
% Display
fprintf('Task 3:\n');
for f_tmp = 1:task3_samples
    fprintf('When f = %f: ', task3_f(f_tmp));
    fprintf('F1 = %f + %fi\n', task3_F1(f_tmp), task3_F1(f_tmp)/1i);
end
disp('---');
    

%% Task 4: Calculate the interpolated results and test their accuracy 
% Local parameters
task4_f = task3_f;
task4_F1 = task3_F1;
task4_f_odd = task3_f(1:2:end);        % samples: retain odd indices frequencies 
task4_F1_odd = task3_F1(1:2:end);      % samples: retain odd indices function values 
% Interpolation
% Note: pchip(x,y,xq) return a vector of interpolated values corresponding to the query points in xq
task4_fldint = pchip(task4_f_odd, task4_F1_odd, task4_f);  
% Calcualte the normalized errors 
% Note: normalized error = abs(interpolated value - function value) / max(abs(max function value))
task4_Error = abs( task4_fldint - task4_F1 ) / max( abs(task4_F1) );
% Display
fprintf('Task 4:\n');
for f_tmp = 1:length(task4_f)
    fprintf('When f = %f:\n', task4_f(f_tmp));
    fprintf('function value F1 = %f + %fi, ', task4_F1(f_tmp), task4_F1(f_tmp)/1i);
    fprintf('interpolated result F1 = %f + %fi, ', task4_fldint(f_tmp), task4_fldint(f_tmp)/1i);
    fprintf('Error = %f\n', task4_Error(f_tmp));
end
disp('Note: Becasue the interpolation must pass through sample points,');
disp('so errors of the functiion values in the odd indices (sample points) are 0.');
disp('---');


%% Task 5: Bisect and obtain new frequency points through loop
% Local parameters
task5_f = task3_f;             % Initial: its items will change (add) with the below loop
task5_Error = task4_Error;     % Initial: its items will change (add) with the below loop
max_error = 0.01;              % the condition: max error 
error_point = [];              % Record the points with an error greater than max_error for each round of the while loop
num = 0;                       % Counter: record the total number of the while loop
% Dual Loop: Outer loop is while-loop, the inner has 3 for-loop. 
% Note1: We only need to check even points, cuz odd points must be known as samples in later interpolation.
% Note2: When we find an even point doesn't meet the condition, 2 new points are added around (left and right) this point,
%        and we still retain this even point, which means: 
%        the indice of the old even point becomes odd! while the indices of the 2 new added points are even!
while max(task5_Error) > max_error
    % First, check all error points in task5_f thoroughly, and record their frequencies into error_point
    for even_tmp = 2:2:length(task5_f)                        % only check even points
        if task5_Error(even_tmp) > max_error                  % meet the 'bad' condition
            error_point = [error_point, task5_f(even_tmp)];   % record its frequency
        end
    end
    % Second, determine if we've found new error points, if so we continue to add new frequencies in task5_f
    if isempty(error_point) ~= 1
        for error_tmp = error_point  
            indice_tmp = find(task5_f == error_tmp);    % the indice of error point in current task5_f
            indice_tmp_left = indice_tmp - 1;           % the indice of the point to the left of the error point
            indice_tmp_right = indice_tmp + 1;          % the indice of the point to the right of the error point
            % Calculate the 2 new frequencies and add them into task5_f
            f_tmp_left = 10^( (log10(task5_f(indice_tmp_left)) + log10(task5_f(indice_tmp))) / 2 );     % new freq1 (left)
            f_tmp_right = 10^( (log10(task5_f(indice_tmp)) + log10(task5_f(indice_tmp_right))) / 2 );   % new freq2 (right)
            % Adding method: f = [odd_left, new1, error_point, new2, odd_right]
            task5_f = [ task5_f(1:indice_tmp_left), f_tmp_left, error_tmp, f_tmp_right, task5_f(indice_tmp_right:end) ];
        end
    else
        % If we've not found new error points, we can end the while-loop directly,
        % which means current all interpolation errors of all frequency points are less than 0.01!
        % so we don't need to add new interpolation points, and we can end the whole loop.
        break   
    end
    % Third, use the new task5_f and calculate new F1 (integration) on these new frequencies. 
    task5_F1 = zeros(1, length(task5_f));
    for f_tmp = 1:length(task5_f) 
        omega_tmp = 2 * pi * task5_f(f_tmp);   
        gamma_tmp = sqrt( 1i * omega_tmp * sigma * mu );
        Integrand_tmp = @(kappa) Integrand_Func(kappa, gamma_tmp, rho, h);  
        F1_tmp = 1i * omega_tmp * quadgk(Integrand_tmp, 0, inf); 
        task5_F1(f_tmp) = F1_tmp;
    end
    % Fourth, select odd points from the new task5_f as new samples (interpolation) and calculate the new errors. 
    task5_f_odd = task5_f(1:2:end);
    task5_F1_odd = task5_F1(1:2:end); 
    task5_fldint = pchip(task5_f_odd, task5_F1_odd, task5_f);
    task5_Error = abs( task5_fldint - task5_F1 ) / max( abs(task5_F1) );
    % Fifth, we must clean the error_point after a cycle, cuz the error points must be different in next cycle.
    error_point = [];
    num = num + 1;    % counter + 1
end
% Plot
task5_re_part = real(task5_fldint);       % the real part of interpolated values
task5_im_part = imag(task5_fldint);       % the imaginary part of interpolated values
Figure1 = figure();                       % open a new figure 
semilogx(task5_f, task5_re_part, 'x');    % semilogx sets the x-axis use a base-10 logarithmic scale
hold on;                                  % plot the real and imaginary parts into one figure 
semilogx(task5_f, task5_im_part, 'o');
legend('real','imaginary')
title('A3.2-T5 Real and Imaginary Parts of the Interpolated Values');
xlabel('Frequency (Hz)');
ylabel('Space Frequency Domain Function');
% Display
disp('Task 5:');
fprintf('After %d loops of bisecting, we finally need %d frequencies to compute the integral.\n', num, length(task5_f));
disp('In figure 1, we can obviously see the x-axis (frequency axis) is discretized and irregular.');
disp('---');


%% Task 6: Compute the field for all linear frequencies 
%          by interpolating the field computed at the limitd irregular frequencies obtained in Task 5
% Local parameters
task6_f_start = f_start;      % the start frequency is 10^(-2.5) Hz
task6_f_end = f_end;          % the end frequency is 10^(2.5) Hz
task6_f_space = f_start;      % space = f_start = 10^(-2.5)
nf = (task6_f_end - task6_f_start) / task6_f_space + 1;   % the total number of frequencies
task6_f = task6_f_start : task6_f_space : task6_f_end;
% Interpolate the field with new interval (task6_f)
task6_fldint = pchip(task5_f, task5_fldint, task6_f);
task6_re_part = real(task6_fldint);       % the real part of interpolated values
task6_im_part = imag(task6_fldint);       % the imaginary part of interpolated values
% Plot
Figure2 = figure();                       % open a new figure 
semilogx(task6_f, task6_re_part, 'x');    % semilogx sets the x-axis use a base-10 logarithmic scale
hold on;                                  % plot the real and imaginary parts into one figure 
semilogx(task6_f, task6_im_part, 'o');
legend('real','imaginary')
title('A3.2-T6 Real and Imaginary Parts of the Interpolated Values');
xlabel('Frequency (Hz)');
ylabel('Space Frequency Domain Function');


%% Task 7:
disp('Task 6:');
disp("We don't fully understand the task 7. We have no idea about this question.");
disp('---');
dt = 1 / (task6_f_space*2*nf);
task7_t = zeros(1, int64(nf));
task7_t(1) = 0.01;
for t_tmp = 1:nf
    task7_t(t_tmp+1) = task7_t(1) + dt * (t_tmp-1);
end
% IFT
awint = task6_fldint;
awt = 2 * real( ifft(awint) ) / dt;
% Plot
Figure3 = figure();
% x7_range = logspace( log10(task7_t(1)), log10(task7_t(end)), nf+1);
x7_range = task7_t;
y7_range = awt;
loglog(x7_range, y7_range);






