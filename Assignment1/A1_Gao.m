%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: AESM1511-Assignment1 
% Name: Boyu Gao
% Student Number: 5794676
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all parameters 
clear all; clc;    

% 6 scales
scale = [1, 3, 5, 7, 9, 11];

% horizontal axis (x-axis): sampling points 
total_points = 1257;
start_value = -pi/2;
end_value = 3*pi/2;
x = linspace(start_value, end_value, total_points);
fprintf('Total points along the horizontal axis: %d\n', total_points);
fprintf('Sampling interval: %f\n', (end_value - start_value)/(total_points-1) );

% vertical axis (y-axis): 6 signals 
y1 = sin( x*scale(1) ) / scale(1);  % the master signal
y2 = sin( x*scale(2) ) / scale(2);  % signal 2 
y3 = sin( x*scale(3) ) / scale(3);  % signal 3
y4 = sin( x*scale(4) ) / scale(4);  % signal 4
y5 = sin( x*scale(5) ) / scale(5);  % signal 5
y6 = sin( x*scale(6) ) / scale(6);  % signal 6
y = [y1; y2; y3; y4; y5; y6];       % generate a matrix: 6 x 1257; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task4 - the sum of the first 2 / 3 / 4 / 5 / 6 signals  
master = y(1,:);                                                % master 
first2 = y(1,:) + y(2,:);                                       % the sum of the first 2
first3 = y(1,:) + y(2,:) + y(3,:);                              % the sum of the first 3
first4 = y(1,:) + y(2,:) + y(3,:) + y(4,:);                     % the sum of the first 4
first5 = y(1,:) + y(2,:) + y(3,:) + y(4,:) + y(5,:);            % the sum of the first 5
first6 = y(1,:) + y(2,:) + y(3,:) + y(4,:) + y(5,:) + y(6,:);   % the sum of the first 6
% figure1: plot
task4_fig = figure(1);
plot(x, master, 'DisplayName', 'master'); hold on;  
plot(x, first2, 'DisplayName', 'first2'); hold on;  
plot(x, first3, 'DisplayName', 'first3'); hold on;  
plot(x, first4, 'DisplayName', 'first4'); hold on;  
plot(x, first5, 'DisplayName', 'first5'); hold on;  
plot(x, first6, 'DisplayName', 'first6'); 
% add legend
lgd = legend;
% other set
grid on; xlabel('x (radian)');  ylabel('y'); title('Figure1-Task4');
xlim([start_value, end_value]);  % restrict the x-axis range
% save figure (.fig)
figure(task4_fig)
savefig('Task4.fig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task5 - the sum of even-number / odd-number / even-number + odd_number matrix rows  
y_odd = y(1,:) + y(3,:) + y(5,:);   % the sum of odd-number matrixs rows 
y_even = y(2,:) + y(4,:) + y(6,:);  % the sum of even_number matrix rows 
y_odd_even = y_odd + y_even;        % the sum of odd-number + even_number matrix rows 
% figure2: plot
task5_fig = figure(2);
plot(x, y_odd, 'DisplayName', 'sum of odd'); hold on;
plot(x, y_even, 'DisplayName', 'sum of even'); hold on;
plot(x, y_odd_even, 'DisplayName', 'sum of odd + even');
% legend
lgd = legend;
% other set
grid on; xlabel('x (radian)');  ylabel('y'); title('Figure2-Task5');
xlim([start_value, end_value]);  % restrict the x-axis range
% save figure (.fig)
figure(task5_fig)
savefig('Task5.fig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure3: task6 - two finals summed results of task4 and task5
task4 = first6;
task5 = y_odd_even;
% figure3: plot
task6_fig = figure(3);
plot(x, task4, 'DisplayName', 'task4'); hold on;
plot(x, task5, 'DisplayName', 'task5');
% legend
lgd = legend;
% other set
grid on; xlabel('x (radian)'); ylabel('y'); title('Figure3-Task6');
xlim([start_value, end_value]);  % restrict the x-axis arange
% save figure (.fig)
figure(task6_fig)
savefig('Task6.fig');