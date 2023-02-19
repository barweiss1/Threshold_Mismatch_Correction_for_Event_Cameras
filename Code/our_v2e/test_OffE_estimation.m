%%  A test for the OffE calibration method
%   Asaf Omer and Bar Weiss
% 
% The code below implements a v2e convertion and the OffE threshold
% estimation method on resulting stream.
%
% The OffE method is taken from Ziwei Wang et al "Event Camera Calibration
% of Per-Pixel Biased Contrast Threshold".
% 
%

%%  Setup

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
fontSize = 14;
   
% Change the current folder to the folder of this m-file.
% (The line of code below is from Brett Shoelson of The Mathworks.)
workingDir = fileparts(matlab.desktop.editor.getActiveFilename);
cd(workingDir);
workingDir_str = convertCharsToStrings(workingDir);
imagesDir_str = workingDir_str + '\images';

rmdir(imagesDir_str ,'s')
mkdir('images')

addpath('../threshold correction');
  
%video_to_use_str = 'generated videos/our_video.mp4'; %set name of video
video_to_use_str = 'generated videos/flicker_video.avi'; %set name of video
%video_to_use_str = 'generated videos/bars_vertical_video.avi'; %set name of video
%video_to_use_str = 'generated videos/bars_horizontal_video.avi'; %set name of video
%video_to_use_str = 'generated videos/moving_sine_video.avi'; %set name of video
%video_to_use_str = 'generated videos/generated videos/linear_gradient_video.avi'; %set name of video
%video_to_use_str = 'generated videos/linear_gradient_video_8s.avi'; %set name of video

movieFullFileName = fullfile(workingDir, video_to_use_str);
video = VideoReader(movieFullFileName,'CurrentTime',0);

%% Generate frames from given video

frames_no = 0;
frames_array = cell([],1) ;
gs_frames_array = cell([],1);
mat_frames_array = cell([],1);

while hasFrame(video)
    curr_frame = readFrame(video);    
    frames_no = frames_no + 1;
    frames_array{frames_no} = curr_frame ;
    gs_frames_array{frames_no} = im2double(rgb2gray(curr_frame))*255;
end

frame_size = size(curr_frame);
frame_size_x = frame_size(2);
frame_size_y = frame_size(1);

%% Below are several threshold patterns that can be selected to understand the threshold mismatch effect
%% uniform threshold matrix

positive_threshold_matrix = ones([frame_size_y,frame_size_x])*5;
negative_threshold_matrix = ones([frame_size_y,frame_size_x])*-5;

%% asymmetric - normal distribution threshold matrix

pos_mu = 0.3;
neg_mu = -0.3;

pos_sigma = 0.03;
neg_sigma = 0.03;

positive_threshold_matrix = normrnd(pos_mu,pos_sigma,[frame_size_y,frame_size_x]);
negative_threshold_matrix = normrnd(neg_mu,neg_sigma,[frame_size_y,frame_size_x]);

%% symmetric - normal distribution threshold matrix

% the negative matrix needs to be negative 

pos_mu = 0.15;
neg_mu = -0.15;

pos_sigma = 0.03;
neg_sigma = 0.03;

positive_threshold_matrix = normrnd(pos_mu,pos_sigma,[frame_size_y,frame_size_x]);
negative_threshold_matrix = -1*positive_threshold_matrix;

% % first 20 raws as reference
% positive_threshold_matrix(1:20,:) = pos_mu;
% negative_threshold_matrix(1:20,:) = neg_mu;



%% Create Events frames by subtracting subsequent frames

rgbImage_event = cell([],1);

events_timestamp = [];
events_x = [];
events_y = [];
events_polarity = [];

output_fps = 25;
input_fps = video.FrameRate;
fps_ratio = floor(input_fps/output_fps);
fps_counter = 0;

%pixel_values = zeros(size(gs_frames_array{1}));
pixel_values = linlog(gs_frames_array{1});

pos_events_video = zeros(size(gs_frames_array{1}));
neg_events_video = zeros(size(gs_frames_array{1}));

for frame = 1:(frames_no - 1)
    
   % calculate log-ilumination changes
   %diff_frames = log(gs_frames_array{frame + 1}+0.001) - log(gs_frames_array{frame}+0.001); 
   diff_frames = linlog(gs_frames_array{frame + 1}) - linlog(gs_frames_array{frame}); 
   pixel_values = pixel_values + diff_frames;
   
   % if the change is higher than the threshold than generate an event 
   positive_events_map = pixel_values > positive_threshold_matrix;
   negative_events_map = pixel_values < negative_threshold_matrix;
   
   % count the number of events at each pixel
   positive_events = floor(positive_events_map.*pixel_values./positive_threshold_matrix);
   negative_events = floor(negative_events_map.*pixel_values./negative_threshold_matrix);
   
   % reset intensity values for triggered events
   pixel_values = pixel_values - positive_events.*positive_threshold_matrix;
   pixel_values = pixel_values - negative_events.*negative_threshold_matrix;
   
   total_events = positive_events + negative_events;
   
   % save timestamps 
   events_num = sum(sum(positive_events)) + sum(sum(negative_events));
   events_timestamp = [events_timestamp repmat((frame/video.FrameRate)*1e6,[1 events_num])];
   % seperate positive and negative event indecies so in order to save
   % polarities in order
   [y_pos_idx, x_pos_idx] = map2xy_stream(positive_events);
   [y_neg_idx, x_neg_idx] = map2xy_stream(negative_events);
   
   events_x = [events_x (x_pos_idx' - 1) (x_neg_idx' - 1)];
   events_y = [events_y (y_pos_idx' - 1) (y_neg_idx' - 1)];
   
   % we entered the positive polarities first
   events_polarity = [events_polarity ones(1,sum(sum(positive_events))) zeros(1,sum(sum(negative_events)))];
   
   % add current events to video interpolation
   pos_events_video = pos_events_video | positive_events_map;
   neg_events_video = neg_events_video | negative_events_map;
   
   % interpolate to output frame rate for event video format 
   fps_counter = fps_counter + 1;
   if fps_counter >= fps_ratio
      
      fps_counter = 0;
      % convert to color values
      positive_events_color = pos_events_video*255;
      negative_events_color = neg_events_video*255;

      % reset interpolation images
      pos_events_video = zeros(size(positive_events_color));
      neg_events_video = zeros(size(negative_events_color));
      
      % add to videos
      rgbImage_event = cat(3, negative_events_color,positive_events_color,positive_events_color*0);

      filename = [sprintf('%03d',frame) '.jpg'];
      fullname = fullfile(workingDir,'images',filename);
      imwrite(rgbImage_event,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
      
   end

   
end



%% Make video out of created frames

outputVideo = VideoWriter(fullfile(workingDir,'generated videos/linear_gradient_events_8s.avi'));
%outputVideo.FrameRate = video.FrameRate;
outputVideo.FrameRate =30;
open(outputVideo)

imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = {imageNames.name}';

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)

%% Single Pixel test

% define start and end time
Tstart = frames_no * 0.2 * 1e6 / video.FrameRate;
Tend = frames_no * 1 * 1e6 / video.FrameRate;
window_num = floor(frames_no*0.8/fps_ratio);

% select pixel coordinates
x_0 = 85;
y_0 = 124;

thresh_pos = positive_threshold_matrix(y_0,x_0);
thresh_neg = negative_threshold_matrix(y_0,x_0);

% seperate one pixel stream
[ts_cell,pol_cell] = stream_map(events_x,events_y,events_timestamp,events_polarity,frame_size_x,frame_size_y);

ts = ts_cell{y_0,x_0};
pol = pol_cell{y_0,x_0};


T = linspace(Tstart,Tend,window_num);

[m,rms] = single_pix_thresh(ts,pol,T,true);

% repeat for 2nd pixel
x_1 = 201;
y_1 = 30;
ts1 = ts_cell{y_1,x_1};
pol1 = pol_cell{y_1,x_1};
thresh_pos1 = positive_threshold_matrix(y_1,x_1);
thresh_neg1 = negative_threshold_matrix(y_1,x_1);

% plot signals
t = linspace(Tstart,Tend,10000);
zoh = ZOH_rec(ts,pol,thresh_pos,thresh_neg,t);
zoh1 = ZOH_rec(ts1,pol1,thresh_pos1,thresh_neg1,t);
figure;
plot(t,zoh);
hold on
plot(t,zoh1);
title('ZOH reconstruction');
legend('1st pixel','2nd pixel');

%% Full Threshold Estimation

[c_est,b_est] = OffE_thresh_est(events_x,events_y,events_timestamp,events_polarity,window_num,Tstart,Tend,frame_size_x,frame_size_y);

% convert to positive and negative threshold format
thresh_pos_est = c_est + b_est;
thresh_neg_est = c_est - b_est;

%% plot

% compare positive thresholds 
% normalize for comparison
norm_est_pos = thresh_pos_est/max(max(thresh_pos_est));
norm_real_pos = positive_threshold_matrix/max(max(positive_threshold_matrix));
diff = (norm_est_pos - norm_real_pos)./norm_real_pos;
figure(4)
imagesc(diff);
colorbar
title('estimation error')
avg_error = sum(sum(abs(diff)))/numel(diff);
display("The positive threshold MAE estimation error is " + avg_error*100 + "%")
avg_error = sqrt(sum(sum(diff.^2))/numel(diff));
display("The positive threshold MSE estimation error is " + avg_error*100 + "%")

figure(3)
subplot(1,2,1);
imagesc(norm_est_pos)
colorbar
title('estimated positive thresholds')
subplot(1,2,2);
imagesc(norm_real_pos)
colorbar
title('actual positive thresholds')

% compare negative thresholds 
% normalize for comparison
norm_est_neg = thresh_neg_est/max(max(abs(thresh_neg_est)));
norm_real_neg = abs(negative_threshold_matrix)/max(max(abs(negative_threshold_matrix)));
diff = (norm_est_neg - norm_real_neg)./norm_real_neg;
avg_error = sum(sum(abs(diff)))/numel(diff);
display("The negative threshold MAE estimation error is " + avg_error*100 + "%")
avg_error = sqrt(sum(sum(diff.^2))/numel(diff));
display("The negative threshold MSE estimation error is " + avg_error*100 + "%")

figure(5)
subplot(1,2,1);
imagesc(norm_est_neg)
colorbar
title('estimated negative thresholds')
subplot(1,2,2);
imagesc(norm_real_neg)
colorbar
title('actual negative thresholds')

figure(6)
subplot(2,2,1);
histogram(norm_est_pos,100);
title('Histogram of Estimated Positive Threshold Values');
xlabel('Estimated Threshold Value');
subplot(2,2,2);
histogram(norm_est_neg,100);
title('Histogram of Estimated Negative Threshold Values');
xlabel('Estimated Threshold Value');
subplot(2,2,3);
histogram(norm_real_pos,100);
title('Histogram of Actual Positive Threshold Values');
xlabel('Actual Threshold Value');
subplot(2,2,4);
histogram(norm_real_neg,100);
title('Histogram of Actual Negative Threshold Values');
xlabel('Actual Threshold Value');


