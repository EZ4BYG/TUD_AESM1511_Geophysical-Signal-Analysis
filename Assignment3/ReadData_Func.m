% Function file: Read original reflection common-source gather data.
% Input: fpath, traces, time_samples
% Return: data_refl_time

function [data_refl_time] = ReadData_Func(fpath, traces, time_samples, dat_format, mac_format)
    fid = fopen(fpath, 'r');    % open the file 
    temp = fread(fid, traces * time_samples, dat_format, 0, mac_format);  % fread read the binary file
    fclose(fid);                % close the file
    data_refl_time(:,:) = reshape(temp, time_samples, traces);  % get the time-domain data and reshape its size to 1001 x 401
end