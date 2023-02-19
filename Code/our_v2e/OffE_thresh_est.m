function [c,b] = OffE_thresh_est(x,y,ts,pol,window_num,Tstart,Tend,pixel_num_x,pixel_num_y)

% This function implements the OffE method threshold estimation in
% from "Event Camera Cabliration of Per-pixel Biased Contrast
% Threshold". This function operates on a single pixel and applies curve
% fitting of the sum of polarities to the total event number and calculate
% the RMS of the estimation. 
%
% Inputs - 
% x = horizontal pixel index vector of occuring events, i.e x(i) is the x
% index of the i-th event
% y = vetical pixel index vector of occuring events, i.e y(i) is the y
% index of the i-th event
% ts = timestamp vector of occuring events 
% pol = polarity vector of events 1 is ON event, 0 is OFF event
% window_num = number of windows the measring time needs to be divided to
% Tstart = start time of the scene 
% Tend = end time of the scene 
% pixel_num_x = number of pixels in x axis
% pixel_num_y = number of pixels in y axis
% 
% Outputs - 
% c = a matrix containing the estimated contrast threshold for each pixel 
% b = a matrix containing the estimated bias of the conrast threshold for
% each pixel 
% 
%

% ----- code -----

% convert streams to per pixel format
[ts_cell, pol_cell] = stream_map(x,y,ts,pol,pixel_num_x,pixel_num_y);

% vector of timeframes
T = linspace(Tstart,Tend,window_num);

% calculate fit and rms for each pixel
wrapper = @(s,r) single_pix_thresh(s,r,T,false); % define function with shared T vector input for cellfun
[m,rms] = cellfun(wrapper,ts_cell,pol_cell,'UniformOutput',true);


% estimate according to the paper
rms_med = median(rms(~isnan(rms))); % median is calculated just for the valid results
c = rms_med./rms;
b = -c.*m;


% replace invalid values with the default threshold 
replace_mask = isnan(rms) | isnan(c) | isnan(b) | (abs(c)>10);
c(replace_mask) = mean(c(~replace_mask));
b(replace_mask) = 0;

end