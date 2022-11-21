%This is an example solution to Assignment 3.1
clear variables
close all
clc 
%% Task 1
%Input parameters of the reflection data
filename_refl_data='refl_3layers_fp50_dx0p5_500_rvz_tapered.bin'; %file with reflection common-source gather
read_prec='float32'; %precision for reading the binary values from the file
read_form='ieee-le'; %format for reading the binary values from the file
traces=401; %the number of traces in the input file
ns=1001; %the number of time samples in each trace
dt=0.001; %the time sampling in each trace in seconds
dx=2.5; %the spatial sampling in meters, i.e., the distance between two neighboring traces
xfirsttrace=510; %the position of the first receiver along the line

fprintf(['The time sampling is ',num2str(dt,'%5.4f'),' (s). \n'])
fprintf(['The spatial sampling is ',num2str(dx,'%5.1f'),' (m). \n'])
fprintf(['The numer of traces is ',num2str(traces,'%5.0f'),'. \n'])

%Parameters used in the frequncy-space domain
fcut=200; %frequencies higher than this one are expected to be zero and are discarded

%Parameters used in the frequncy-wavenumber domain
ntaper=10; %points to taper along the wavenumber axis before setting values to zero
taperexcur=20; %points to the left (or right) of the surface-wave energy maximum to start taper
v_sw=450; %the estimated highest surface-wave velocity

%Plotting parameters
fs=14; %Fontsize
lw=2; %Linewidth
cmap=gray(256); %the colour map for plotting the common-source gather
dispscale=0.01; %0.01; %scaling the colorbar to clip the data to display weaker arrivals
dispscalefilt=0.01; %0.1; %scaling the colorbar to clip the data to display the filtered result

%Allocating memory for the reflection data
fprintf('Allocate memory for the data matrix...')
data_refl=zeros(ns,traces);
fprintf('done\n')

%Now we are going to read the binary file into Matlab using the function
data_refl(:,:)=readbinsu(filename_refl_data,read_prec,read_form,traces,ns);

%Defining the time and space axes
fprintf('Make the axis vectors...')
timeaxis=linspace(0,ns-1,ns)*dt;
spaceaxis=linspace(0,traces-1,traces)*dx+xfirsttrace;
fprintf('done\n')

figure(1);
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
%Transforming the data to the frequency-wavenumber domain

fprintf('To transform the data to the frequency domain, we are using fft() \n')
fprintf('because it is defined as "exp(-i*2*pi*f)". \n')
fprintf('To transform the data to the wavenumber domain, we are using ifft() \n')
fprintf('because it is defined as "exp(i*2*pi*k)". \n')

fprintf('FFT to frequency-space domain...')
data_refl_f=fft(data_refl,[],1)*dt;
df=1/(ns*dt); %the frequency sampling in Hertz
nf=ns; %the number of frequency samples
freqaxis=linspace(0,nf/2,nf/2+1)*df; %making the frequency axis
fprintf('done\n')

fprintf('Extract relevant frequencies...')
fcutel=find(freqaxis>=fcut,1,'first');
freqaxiscut=freqaxis(1:fcutel); %the new frequency axis
data_refl_f=data_refl_f(1:fcutel,:);
fprintf('done\n')

figure(2);
imagesc(spaceaxis,freqaxiscut,abs(data_refl_f(:,:)));
colorbar
%caxis([dispscale*min(min(abs(data_refl_f(:,:)))) dispscale*max(max(abs(data_refl_f(:,:))))])
xlabel('Horizontal Distance (m)','Fontsize',fs)
ylabel('Frequency (Hz)','Fontsize',fs)
title('Reflection data in the frequency domain - Amplitude','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

fprintf('IFFT to frequency-wavenumber domain...')
data_refl_fk=fftshift(traces*ifft(data_refl_f,[],2)*dx,2);
dk=1/(traces*dx); %the wavenumber sampling in 1/m
nk=traces; %the number of wavenumber samples
kaxis=linspace(-nk/2,nk/2-1,nk)*dk;
fprintf('done\n')

figure(3);
imagesc(kaxis,freqaxiscut,abs(data_refl_fk(:,:)))
colorbar
xlabel('Wavenumber (1/m)','Fontsize',fs)
ylabel('Frequency (Hz)','Fontsize',fs)
title('Reflection data in the frequency-wavenumber domain','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

%% Task 3
%Suppressing (tapering) the surface waves

%looping over frequencies to find the right wavenumber to taper to the left
%or right of it
data_refl_fkfil=data_refl_fk;
%checking which is the maximum frequency (element) corresponding to the 
%maximum wavenumber (element) present
fswel_max=find(freqaxis<=(kaxis(nk)*v_sw),1,'last');
%looping only till the frequency that ensures we have a corresponding 
%pair (frequency,wavenumber)
for freq=1:(min(fswel_max,fcutel)-1)
    %below we calculate a wavenumber value using the current frequency and the 
    %given surface-wave velocity
    ksw=freqaxiscut(freq)/v_sw;
    kswel_pos=find(kaxis>=ksw,1,'first');
    kswel_neg=find(kaxis<=-ksw,1,'last');
    %tapering values along the wavenumber axis starting from kswel_pos 
    %shifted to the left by taperexcur and going to the right
    lengthk=nk-kswel_pos+taperexcur;
    lintaper=[linspace(1,0,ntaper) zeros(1,lengthk-ntaper)];
    data_refl_fkfil(freq,(kswel_pos-taperexcur+1):end)=data_refl_fk(freq,(kswel_pos-taperexcur+1):end).*lintaper;
    %tapering values along the wavenumber axis starting from kswel_neg 
    %shifted to the right by taperexcur and going to the left
    lengthk=kswel_neg+taperexcur;
    lintaper=[zeros(1,lengthk-ntaper) linspace(0,1,ntaper)];
    data_refl_fkfil(freq,1:(kswel_neg+taperexcur))=data_refl_fk(freq,1:(kswel_neg+taperexcur)).*lintaper;
    
end

figure(4);
imagesc(kaxis,freqaxiscut,abs(data_refl_fkfil(:,:)))
colorbar
xlabel('Wavenumber (1/m)','Fontsize',fs)
ylabel('Frequency (Hz)','Fontsize',fs)
title('F-k filtered reflection data in the frequency-wavenumber domain','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

%Transforming the data back to the frequency-space domain
fprintf('FFT to frequency-space domain...')
tempfdata=fft(fftshift(data_refl_fkfil,2),[],2)*dk;
fprintf('done\n')

%Adding zeroes in stead of the removed frequencies to obtain time-domain
%data as long as before the transformations.
fprintf('Add removed frequencies...')
addfreq = zeros(nf,traces);
addfreq(1:fcutel,:) = tempfdata;
fprintf('done\n')

%Transforming the data back to the time-space domain
fprintf('IFFT to time-space domain...')
data_refl_swfilt=2*real(nf*ifft(addfreq,[],1)*df);
fprintf('done\n')

figure(5);
imagesc(spaceaxis,timeaxis,(data_refl_swfilt(:,:)));
colormap(cmap)
colorbar
caxis([dispscalefilt*min(min(data_refl_swfilt(:,:))) dispscalefilt*max(max(data_refl_swfilt(:,:)))])
xlabel('Horizontal Distance (m)','Fontsize',fs)
ylabel('Two-way traveltime (s)','Fontsize',fs)
title('F-k filtered reflection data','Fontsize',fs)
set(gca,'Fontsize',fs)
set(gca,'LineWidth',lw)

fprintf('We can see that the f-k filtering has suppressed the surface-wave \n')
fprintf('energy, while the body-wave energy has been preserved. \n')
fprintf('We can also notice that there are artefacts from the filtering. \n')


