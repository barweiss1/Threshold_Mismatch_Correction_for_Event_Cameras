function [freqs,threshold] = thresh_1_pix(x,y,ts,pol,ts_ref,pol_ref,light_function)

% This function estimates the threshold of each pixel. this is done using
% the relation - Event rate(t) = (1/threshold)*(dln(I)/dt). if we know the
% light pattern than we can use that to extract the threshold
%
% Inputs - 
% x = horizontal pixel index vector of occuring events, i.e x(i) is the x
% index of the i-th event
% y = vetical pixel index vector of occuring events, i.e y(i) is the y
% index of the i-th event
% ts = timestamp vector of occuring events , i.e ts(i) is the timestamp of
% the i-th event
% pol = polarity vector of events 1 is ON event, 0 is OFF event
% pixel_num_x = number of pixels in x axis
% pixel_num_y = number of pixels in y axis
% light_function = a handle (pointer) to the fucntion that represents the lighting
% change in time
% 
% Outputs - 
% threshold = a matrix with the estimated threshold of each pixel
% 
%

% define cells to save running frequency estimation
freqs = [];
threshold = [];

prev_ts = ts(1);

%a = light_function(ts(1),x,y);
%b = light_function(ts(2)/2,x,y);
%c = light_function(ts(3),x,y);

for i=2:length(ts)
    
    % estimate event frequency
    f_est = 1/(ts(i) - prev_ts);
    freqs = [freqs f_est];
    
    % get the threshold from the lighting function 
    % need to modify this equation to include polarity and diffrent ON and
    % OFF thresholds 
    th_est = abs(log(light_function(ts(i),x,y)) - log(light_function(prev_ts,x,y)));
    threshold = [threshold th_est];

    prev_ts = ts(i);
    
end

end