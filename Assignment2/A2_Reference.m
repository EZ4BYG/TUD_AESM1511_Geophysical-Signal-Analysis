%This is an example solution to Assignment 2
clear variables
close all
clc 
%% Task 1
%%You can put here as general comments some description from a task (in
%%this case Task 1) to give a context to the part that follows. I do not do
%%it as you know what Task 1 contains (i.e., I know the expected audience)

%Input parameters to read the reflection data
filename_refl_data='./refl_3layers_fp50_dx0p5_500_rvz.bin'; %file with reflection common-source gather
read_prec='float32'; %precision for reading the binary values from the file
read_form='ieee-le'; %format for reading the binary values from the file
traces=401; %the number of traces in the input file
ns=1001; %the number of time samples in each trace
dt=0.001; %the time sampling in seconds in each trace
dx=2.5; %the spatial sampling in meters, i.e., the distance between two neighboring traces
xshot=500; %the position of the shot along the line
xfirsttrace=510; %the position of the first receiver along the line

fprintf(['The time sampling is ',num2str(dt,'%5.4f'),' (s). \n'])
fprintf(['The spatial sampling is ',num2str(dx,'%5.1f'),' (m). \n'])
fprintf(['The numer of traces is ',num2str(traces,'%5.0f'),'. \n'])

%Parameters used in the frequncy-space domain
fcut=200; %frequencies higher than this number are expected to be zero and are discarded
flowcut=60; %frequencies below this number will be muted for frequency-filtering purposes
ftaper=30; %frequency above which linear and sin tapers are applied
surf_wave_dist=950; %surface waves present only up to this offset

% Plotting parameters
fs=14; %Fontsize
lw=2; %Linewidth
cmap=gray(256); %the colour map for plotting the common-source gather
dispscale=0.1; %scaling the colorbar to clip the data to display weaker arrivals

% Parameters used to sample the spatial dimention more coarsly for Task 6
xsubsamp=[0 1 2 3]; %this means that every 2^0, 2^1, 2^2, 2^3 trace will be taken

%Allocating memory for the reflection data
fprintf('Allocate memory for the data matrix...')
data_refl=zeros(ns,traces);
fprintf('done\n')

%Now we are going to read the binary file generated under Seismic Unix into Matlab.
fid=fopen(filename_refl_data,'r');
temp=fread(fid,traces*ns,read_prec,0,read_form);
fclose(fid);
data_refl(:,:)=reshape(temp,ns,traces);

%Defining the time and space axes
fprintf('Make the axis vectors...')
timeaxis=linspace(0,ns-1,ns)*dt;
spaceaxis=linspace(0,traces-1,traces)*dx+xfirsttrace;
fprintf('done\n')

figure;
imagesc(spaceaxis,timeaxis,data_refl(:,:));
colormap(cmap)
colorbar
caxis([dispscale*min(min(data_refl(:,:))) dispscale*max(max(data_refl(:,:)))])
xlabel('Horizontal Distance (m)','Fontsize',fs)
ylabel('Two-way traveltime (s)','Fontsize',fs)
title('Recorded reflection data','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

%% Task 2
%Transforming the data to the frequency domain
fprintf('FFT to frequency-space domain...')
data_refl_f=fft(data_refl,[],1)*dt;
df=1/(ns*dt); %the frequency sampling in Hertz
nf=ns; %the number of frequency samples
freqaxis=linspace(0,nf/2-1,nf/2)*df; %making the frequency axis
fprintf('done\n')

fprintf('Extract relevant frequencies...')
fcutel=find(freqaxis>=fcut,1,'first');
freqaxiscut=freqaxis(1:fcutel); %the new frequency axis
data_refl_f=data_refl_f(1:fcutel,:);
fprintf('done\n')

fprintf('Attempting to plot directly the common-source gather in the \n')
fprintf('frequency-space domain gives an error because the data after fft \n')
fprintf('is complex, while imagesc requires real data. \n')

figure;
imagesc(spaceaxis,freqaxiscut,abs(data_refl_f(:,:)));
colorbar
xlabel('Horizontal Distance (m)','Fontsize',fs)
ylabel('Frequency (Hz)','Fontsize',fs)
title('Reflection data in the frequency domain - Amplitude','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)


%% Task 3
fprintf('The dominant energy in the frequency domain corresponds to \n')
fprintf('surface waves in the time domain. \n')

%% Task 4
%Filtering the surface waves

fprintf('Perform filtering: \n')
ix_surf_wave=1+(surf_wave_dist-xfirsttrace)/dx; %the horizontal distance index below which we want to filter surface waves

fprintf(['The surface waves are present only till offset ',num2str(surf_wave_dist,'%5.4f'),' (m), \n'])
fprintf('we do not need to apply the filter for longer offsets. \n')


%builtding the linear and sin low-cut filters
ntaper=flowcut-ftaper+1; %number of points of be tapered
lintaper=[zeros(1,ftaper-1) linspace(0,1,ntaper)]; %defining the linear taper
sintaper=sin(lintaper*pi/2); %defining the sin taper

%creating the matrices that will contain the filtered surface-wave frequencies
data_refl_f_lowcut0=data_refl_f;
data_refl_f_lowcutlin=data_refl_f;
data_refl_f_lowcutsin=data_refl_f;
for x=1:ix_surf_wave
   data_refl_f_lowcut0(1:flowcut,x)=0;
   data_refl_f_lowcutlin(1:flowcut,x)=data_refl_f_lowcutlin(1:flowcut,x).*lintaper';
   data_refl_f_lowcutsin(1:flowcut,x)=data_refl_f_lowcutsin(1:flowcut,x).*sintaper';
end
fprintf('The filtering is done. \n')

% figure;
% imagesc(spaceaxis,freqaxiscut,abs(data_refl_f_lowcut0(:,:)));
% colorbar
% xlabel('Horizontal Distance (m)','Fontsize',fs)
% ylabel('Frequency (Hz)','Fontsize',fs)
% title('Reflection data in the frequency domain - Zero filtering','Fontsize',fs)
% set(gca,'Fontsize',fs)
% set(gca,'LineWidth',lw)
% 
% figure;
% imagesc(spaceaxis,freqaxiscut,abs(data_refl_f_lowcutlin(:,:)));
% colorbar
% xlabel('Horizontal Distance (m)','Fontsize',fs)
% ylabel('Frequency (Hz)','Fontsize',fs)
% title('Reflection data in the frequency domain - Linear filtering','Fontsize',fs)
% set(gca,'Fontsize',fs)
% set(gca,'LineWidth',lw)
% 
% figure;
% imagesc(spaceaxis,freqaxiscut,abs(data_refl_f_lowcutsin(:,:)));
% colorbar
% xlabel('Horizontal Distance (m)','Fontsize',fs)
% ylabel('Frequency (Hz)','Fontsize',fs)
% title('Reflection data in the frequency domain - sin filtering','Fontsize',fs)
% set(gca,'Fontsize',fs)
% set(gca,'LineWidth',lw)

fprintf('IFFT to time-space domain...')
data_refl_lowcut0=2*real(fcutel*ifft(data_refl_f_lowcut0,[],1)*df);
data_refl_lowcutlin=2*real(fcutel*ifft(data_refl_f_lowcutlin,[],1)*df);
data_refl_lowcutsin=2*real(fcutel*ifft(data_refl_f_lowcutsin,[],1)*df);
fprintf('done\n')

figure;
imagesc(spaceaxis,timeaxis,data_refl_lowcut0(:,:));
colormap(cmap)
colorbar
caxis([dispscale*min(min(data_refl(:,:))) dispscale*max(max(data_refl(:,:)))])
xlabel('Horizontal Distance (m)','Fontsize',fs)
ylabel('Two-way traveltime (s)','Fontsize',fs)
title('Reflection data filtered using zeros','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

figure;
imagesc(spaceaxis,timeaxis,data_refl_lowcutlin(:,:));
colormap(cmap)
colorbar
caxis([dispscale*min(min(data_refl(:,:))) dispscale*max(max(data_refl(:,:)))])
xlabel('Horizontal Distance (m)','Fontsize',fs)
ylabel('Two-way traveltime (s)','Fontsize',fs)
title('Reflection data filtered using linear taper','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

figure;
imagesc(spaceaxis,timeaxis,data_refl_lowcutsin(:,:));
colormap(cmap)
colorbar
caxis([dispscale*min(min(data_refl(:,:))) dispscale*max(max(data_refl(:,:)))])
xlabel('Horizontal Distance (m)','Fontsize',fs)
ylabel('Two-way traveltime (s)','Fontsize',fs)
title('Reflection data filtered using sin taper','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

fprintf('The surface waves are suppressed, while reflections are \n')
fprintf('relatively preserved. Still, the lower frequencies of the \n')
fprintf('reflections are damaged, in the case of the linear and sin \n')
fprintf('tapers, or lost, in the case of setting frequencies to zero. \n')
fprintf('In the latter case, we also see that strong ringing appears. \n')

%% Tasks 5 and 6 combined 
%Transforming the data to the frequency-wavenumber domain and plotting
figure;
for x=xsubsamp
    subsamp=2^x;
    nsubplot=1+x;
    tempdata=data_refl_f(:,1:subsamp:end);
    data_refl_fk=fftshift(fft(tempdata,[],2)*(dx*subsamp),2);
    dk=1/((traces/subsamp)*dx*subsamp); %the wavenumber sampling in 1/m
    nk=traces/subsamp; %the number of wavenumber samples
    kaxis=linspace(-nk/2,nk/2-1,nk)*dk;

    %As there are four elements in xsubsamp, we will have 4 subplots.
    subplot(2,2,nsubplot)
    imagesc(kaxis,freqaxiscut,abs(data_refl_fk(:,:)))
    colorbar
    caxis([dispscale*min(min(abs(data_refl_fk(:,:)))) dispscale*max(max(abs(data_refl_fk(:,:))))])
    xlabel('Wavenumber (1/m)','Fontsize',fs)
    ylabel('Frequency (Hz)','Fontsize',fs)
    title('Reflection data in the frequency-wavenumber domain','Fontsize',fs)
    set(gca,'Fontsize',fs)
    set(gca,'LineWidth',lw)
end

fprintf('The surface-wave energy is now concentrated along a line. \n')
fprintf('The reflection energy is concentrated around k=0. \n ')




