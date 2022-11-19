%This is an example solution to Assignment 1
%First, we start with a clean sheet: 
clear variables %remove items from workspace, freeing up system memory;
close all %close any open figures;
clc %and clear the command window
%% task 2
%%You can put here as general comments some description from a task (in
%%this case Task 2) to give a context to the part that follows. I do not do
%%it as you know what Task 2 contains (i.e., I know the expected audience).

l_sig=1257; %the length of the signals we are going to use
sigsN=6; %the number of signals we are going to use
s1_beg=-pi/2; %define the beginning of the first signal
s1_end=3*pi/2; %define the end of the first signal
coeff=[1 3 5 7 9 11]; %the coefficients to create the separate signals

sampling=round((s1_end-s1_beg)/l_sig,3); %this is the sampling interval
disp('Answer to the question in Task 2')
disp(['The sampling interval is: ',num2str(sampling,'%5.4f'),' s.'])
disp(' ')

%% Task 3
sigs=zeros(sigsN,l_sig); %define the matrix size and allocate the matrix in the memory
s1=s1_beg:sampling:s1_end; %these are the values of the first signal along the horizontal axis
sigs(1,:)=sin(s1*coeff(1))/coeff(1); %this is the first signal
sigs(2,:)=sin(s1*coeff(2))/coeff(2); %this is the second signal
sigs(3,:)=sin(s1*coeff(3))/coeff(3); %this is the third signal
sigs(4,:)=sin(s1*coeff(4))/coeff(4); %this is the fourth signal
sigs(5,:)=sin(s1*coeff(5))/coeff(5); %this is the fifth signal
sigs(6,:)=sin(s1*coeff(6))/coeff(6); %this is the sixth signal

%% Task 4
xaxis=1:1:l_sig; %defining for the figrues the horizontal axis to be in samples

figure %visualizing the result of summing the signals from 1 to 6
plot(xaxis,sigs(1,:),'r')
hold on

summed=sigs(1,:)+sigs(2,:);
plot(xaxis,summed,'y')

summed=summed+sigs(3,:);
plot(xaxis,summed,'g')

summed=summed+sigs(4,:);
plot(xaxis,summed,'c')

summed=summed+sigs(5,:);
plot(xaxis,summed,'b')

summed=summed+sigs(6,:);
plot(xaxis,summed,'k')
legend('signal 1','signal 1:2','signal 1:3','signal 1:4','signal 1:5','signal 1:6')
hold off

disp('Answer to the question in Task 4')
disp('By adding the six signals one by one, the original master signal')
disp('looses its sinusoidal shape, becomes steeper, flatter at the top,')
disp('and flatter at the bottom.')

disp('If the number of summed signals increases more and more and if')
disp('the signal coefficients increase using the same rule (i.e.,')
disp('"sin(s1*coeff(big))/coeff(big)"), the obtained final signal would look')
disp('like a box (known as boxcar or rectangular signal).')
disp(' ')

%% Task 5
figure %visualizing the result of summing the results from 1,3,5 and from 2,4,6
plot(xaxis,sigs(1,:),'r')
hold on

summed135=sigs(1,:)+sigs(3,:)+sigs(5,:);
plot(xaxis,summed135,'y')

summed246=sigs(2,:)+sigs(4,:)+sigs(6,:);
plot(xaxis,summed246,'g')

summed2=summed135+summed246;
plot(xaxis,summed2,'k')
legend('signal 1','summed signals 1,3,5','summed signals 2,4,6','summed all signals')
hold off

%% Task 6
figure %comparison of the results of the two summation orders
plot(xaxis,summed,'ro',xaxis,summed2,'b-')
legend('summed signals 1,2,3,4,5,6','summed signal 1,3,5+2,4,6')

disp('Answer to the question in Task 6')
disp('Comparing the two final results in Figure 3, we can conclude that')
disp('the summation order DOES NOT matter in obtaining the final result.')
disp('This is a consequence of the superposition principle.')

disp('Comparing the intermediate results in Figures 1 and 2, we can')
disp('conclude that the summation order DOES matter in obtaining')
disp('an intermediate result.')
disp(' ')



