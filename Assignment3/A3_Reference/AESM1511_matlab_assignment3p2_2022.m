%This is an example solution to Assignment 3.2
clear variables
close all
clc 
%% Input parameters
%General input parameters
sigma=1; %electric conductivity
mu=4e-7*pi; %magnetic permeability
f=0.5; %frequency of the input signal
omega=2*pi*f; %angula frequency
grg=1i*omega*sigma*mu; %this is gamma squared
rho=1e3; %horizontal distance
zdist=5e2; %this is the vertical distance between transmitter and receiver
int_beg=0; %integration with equation 1 starts here
int_end=Inf; %integration with equation 1 end here

%Input parameters needed for Task 3
finit=11; %number of initial log frequencies 
fbeg=-2.5; %beginning of the log-frequency interval
fend=2.5; %end of the log-frequency interval

%Plotting parameters
fs=14; %Fontsize
lw=2; %Linewidth

%% Task 2
%Calculating the integral in equation 1 using the function GBkf created for
%Task 1; we integrate over k, so we keep k as a variable
I0K0=1i*omega*quadgk(@(k)GBkf(k,grg,zdist,rho),int_beg,int_end); 

%Calculating the same answer using the zero-order Bessel functions of the 
%first and second kind
IK0=1i*omega*besseli(0,sqrt(grg)*(sqrt(rho^2+zdist^2)-zdist)/2)...
   *besselk(0,sqrt(grg)*(sqrt(rho^2+zdist^2)+zdist)/2);

%Calculating the error between the two results
errgk=abs(I0K0-IK0)/abs(IK0); 
fprintf(['The realtive error between the two methods is ',num2str(errgk,'%14.13f'),'. \n'])
clear errgk I0K0 IK0 %this is needed for the following tasks

%% Task 3
%Now the adaptive interpolation routine starts. This routine minimizes 
%the number of calls to the integration routine

%building the initial frequency axis
lf=linspace(fbeg,fend,finit);
freq=10.^lf;
%determine the number of initial frequencies
nf=length(freq);
%initialize the field array
I0K0=zeros(1,nf);
% compute the field at all initial frequencies
for jf=1:nf
    grg=2i*pi*freq(jf)*sigma*mu; %here, omega is directly exchanged for 2*pi*f
    I0K0(jf)=2i*pi*freq(jf)*quadgk(@(k)GBkf(k,grg,zdist,rho),0,inf);
end
mxf=max(abs(I0K0)); %compute the maximum value of the function to be evaluated

%% Task 4 and 5
%take half of the values and use them for interpolation
IK0=I0K0(1:2:nf);
%create frequency axis for interpolation
wint=2*pi*freq;
%compute the interpolated function
IK0int=pchip(wint(1:2:nf),IK0,wint);
%start the interpolation loop by checking if interpolation is still necessary; the
%array IK0 is used as a dummy with all values from the previous interpolation
%step, the array IK0int contains all the new values
while length(IK0int) ~= length(IK0)
    is=1; %new counter to fill in new frequency values
    IK0=I0K0; %create a new dummy for the next loop
    I0K0(1)=IK0(1); %the first number is retained always because it is computed
    lfi(1)=lf(1); %the new frequency axis is updated for the first value    
    %Below, we do Task 4 and Task 5 together
    for jf=2:2:nf-1
        if abs(IK0int(jf)-IK0(jf))>1e-2*max(abs(IK0)) && abs(IK0(jf))>1e-2*mxf
        %not accurate enough, split interval in three;
            %compute the function left of jf
            lfi(is+1)=(lf(jf-1)+lf(jf))/2; %compute the logarithmic-power value between jf-1 and jf
            snew=1i*2*pi*10^lfi(is+1); %compute the new frequency value between jf-1 and jf
            grg=snew*sigma*mu; %compute grg at the new frequency
            I0K0(is+1)=snew*quadgk(@(k)GBkf(k,grg,zdist,rho),0,inf); %compute the function value at the new frequency
            %the value at jf was already computed and is now stored in the new array
            lfi(is+2)=lf(jf); %update the logarithmic-power value
            I0K0(is+2)=IK0(jf); %shift the known value one index up
            %compute the function at the point right of jf
            lfi(is+3)=(lf(jf+1)+lf(jf))/2; %compute the new logarithmic-power value between jf and jf+1
            snew=2i*pi*10^lfi(is+3); %compute the new frequency value between jf and jf+1
            grg=snew*sigma*mu; %compute grg at new frequency
            I0K0(is+3)=snew*quadgk(@(k)GBkf(k,grg,zdist,rho),0,inf); %compute the function at the new frequency
            is=is+4; %update the counter
            %the value at jf+1 is going to be skipped, because it was already computed and below is being stored in the new array
            lfi(is)=lf(jf+1); %update the logarithmic-power value
            I0K0(is)=IK0(jf+1); %shift the known value one index up
        else
        %step to the next interval and keep the two adjacent precomputed values
            lfi(is+1:is+2)=lf(jf:jf+1); %update logarithmic-power values
            I0K0(is+1:is+2)=IK0(jf:jf+1); %shift the known values one index up
            is=is+2; %update the counter
        end
    end
    %now we estimate the function at the new frequencies with interpolation of
    %the function at the previously chosen frequencies
    IK0int=pchip(wint,IK0,2*pi*10.^lfi);
    %reset the parameters
    lf=lfi;
    wint=2*pi*10.^lf;
    nf=length(lf);
%return to accuracy check and proceed
end
%The plotting below is the end of Task 5
%now we have accurate data and can plot the result
figure(1);
semilogx(10.^lf,real(I0K0),'ro',10.^lf,imag(I0K0),'bx')
axis tight
set(gca,'fontsize',fs)
legend('real part','imaginary part')
xlabel('frequency [Hz]')
ylabel('i\omega I_0(\xi^-)K_0(\xi^+)')
title(['Accurate discrete function at ',num2str(nf,2),' frequencies.'])

%% Task 6
%now interpolate to a linear frequency axis
IK0=[0 IK0];
wint=[0 wint];
df=wint(2)/(2*pi); %calculating the linear-frequency step using the second element in wint
nf=round(max(freq/df)); %calculating the number of liner frequencies
ff=2*pi*(0:nf-1)*df; %the new angular-frequency values 
IK0int=pchip(wint,IK0,ff);
% % If you want to test the accuracy of the overall result, here are the
% % script lines for determining the exact result.
% % Use ctrl+t on the ten lines below this one to turn them into executable lines.
% grg=2i*pi*ff*sigma*mu;
% aw=2i*pi*ff.*besseli(0,sqrt(grg)*(sqrt(rho^2+zdist^2)-zdist)/2)...
%  .*besselk(0,sqrt(grg)*(sqrt(rho^2+zdist^2)+zdist)/2);
% aw(1)=0;
% figure(3); semilogx(ff,abs(aw))
% axis tight
% set(gca,'fontsize',fs)
% xlabel('frequency [Hz]')
% ylabel('i\omega I_0(\xi^-)K_0(\xi^+)')
% title(['exact function at ',num2str(nf,2),' frequencies'])

%% Task 7
% compute time step
dt=pi/ff(nf);
%compute the time-domain function using inverse Fourier transformation 'ifft'
awt=2*real(ifft(IK0int,2*nf))/dt;
%preparing the plots.
t=(2:nf-1)*dt; % compute positive time-axis for values 0.01 seconds and larger
awt=awt(2:nf-1); % keep only positive times 

figure(2) % open figure(2)
tpl=find(awt>=0); % find for which indices the field is positive
tmn=find(awt<0); % find for which indices the field is negative
awtp=nan(1,nf-2); % fill the positive data array with NaN of proper length
awtm=nan(1,nf-2); % fill the negative data array with NaN of proper length
awtp(tpl)=awt(tpl); % fill the positive array with the correct data points
awtm(tmn)=-awt(tmn); % fill the negative array with the correct data points
% plot the result on a log-log plot with a different color for negative
% and poitive function values
loglog(t,awtp,'r',t,awtm,'b')
axis([t(1) t(end) min(abs(awt)) max(abs(awt))])
legend('positive values','negative values')
xlabel('Time (s)')
ylabel('Field strength')
%end of program