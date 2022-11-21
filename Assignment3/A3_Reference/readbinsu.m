function  readbin = readbinsu(fname,r_prec,r_form,ntr,nt)
%readbinsu: this function reads a binary file that has been created from 
%an su file by stripping the headers
%fname - the name of the bin file
%ntr - number of traces in the bin file
%nt - number of time samples in the bin file
%

%Reading the binary file into Matlab.
fid=fopen(fname,'r');
temp=fread(fid,ntr*nt,r_prec,0,r_form);
fclose(fid);
readbin=reshape(temp,nt,ntr);
end