function [ts_cell,pol_cell] = stream_map(x,y,ts,pol,pixel_num_x,pixel_num_y)

% This function converts the single stream AER format to a 2 cell matricies
% with polarity and timestamp streams for each pixel in the camera
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
% 
% Outputs - 
% ts_cell = a matrix with timestamp streams of each pixel
% pol_cell = a matrix with polarity streams of each pixel
% 
%

% define cells to save running frequency estimation
ts_cell = cell(pixel_num_y,pixel_num_x);
pol_cell = cell(pixel_num_y,pixel_num_x);


for i=1:length(x)
   
    ts_cell{y(i)+1,x(i)+1} = [ts_cell{y(i)+1,x(i)+1} ts(i)];
    pol_cell{y(i)+1,x(i)+1} = [pol_cell{y(i)+1,x(i)+1} pol(i)];
    
end


end