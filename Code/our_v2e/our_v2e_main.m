%%  custom v2e tool
%   Asaf Omer and Bar Weiss
% 
% The code below implements a basic tool that converts a syntheticly generated video stream to
% an event stream.
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
  
%video_to_use_str = 'generated videos/our_video.mp4'; %set name of video
video_to_use_str = 'generated videos/flicker_video.avi'; %set name of video
%video_to_use_str = 'generated videos/bars_vertical_video.avi'; %set name of video
%video_to_use_str = 'generated videos/bars_horizontal_video.avi'; %set name of video
%video_to_use_str = 'generated videos/moving_sine_video.avi'; %set name of video
%video_to_use_str = 'generated videos/linear_gradient_video.avi'; %set name of video

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

%% normal distribution threshold matrix

pos_mu = 0.1;
neg_mu = -0.1;

pos_sigma = 0.01;
neg_sigma = 0.01;

positive_threshold_matrix = normrnd(pos_mu,pos_sigma,[frame_size_y,frame_size_x]);
negative_threshold_matrix = normrnd(neg_mu,neg_sigma,[frame_size_y,frame_size_x]);

%% uniform - normal distribution threshold matrix

% the negative matrix needs to be negative 

pos_mu = 0.25;
neg_mu = -0.25;

pos_sigma = 0.05;
neg_sigma = 0.05;

positive_threshold_matrix = normrnd(pos_mu,pos_sigma,[frame_size_y,frame_size_x]);
negative_threshold_matrix = -1*positive_threshold_matrix;

% % first 20 raws as reference
% positive_threshold_matrix(1:20,:) = pos_mu;
% negative_threshold_matrix(1:20,:) = neg_mu;


%% Half and Half

pos_mu_left = 5;
neg_mu_left = -5;
pos_mu_right = 10;
neg_mu_right = -10;

pos_sigma_left = 0.1;
neg_sigma_left = 0.1;
pos_sigma_right = 0.1;
neg_sigma_right = 0.1;

positive_threshold_matrix_left = normrnd(pos_mu_left,pos_sigma_left,[frame_size_y,frame_size_x/2]);
negative_threshold_matrix_left = normrnd(neg_mu_left,neg_sigma_left,[frame_size_y,frame_size_x/2]);
positive_threshold_matrix_right = normrnd(pos_mu_right,pos_sigma_right,[frame_size_y,frame_size_x/2]);
negative_threshold_matrix_right = normrnd(neg_mu_right,neg_sigma_right,[frame_size_y,frame_size_x/2]);

positive_threshold_matrix = [positive_threshold_matrix_left positive_threshold_matrix_right];
negative_threshold_matrix = [negative_threshold_matrix_left negative_threshold_matrix_right];

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

outputVideo = VideoWriter(fullfile(workingDir,'generated videos/flicker_events_1k_fps.avi'));
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


%% Using Histogram / heat map functions
Tstart = frames_no * 0.75 * 1e6 / video.FrameRate;
Tend = frames_no * 0.8 * 1e6 / video.FrameRate;
Xstart = floor(frame_size_x * 0.7);
Xend = floor(frame_size_x * 0.8);
Ystart = floor(frame_size_y * 0.2);
Yend = floor(frame_size_y * 0.9);
nbins = 30;
Fignum = 1;

[hist, heatmap] = timeframe_histogram(events_x,events_y,events_timestamp,events_polarity,Tstart,Tend,Xstart,Xend,Ystart,Yend,nbins,Fignum);
figure(2)
imagesc(heatmap);
colorbar;

%% Threshold Estimation 1 pixel

chosen_pixel_x = frame_size_x/2;
chosen_pixel_y = frame_size_y/2;

chosen_pixel_events_timestamps = [];
chosen_pixel_events_polarity = [];

ref_pixel_events_timestamps = [];
ref_pixel_events_polarity = [];

light_function = @t_linear;
%t_axis = linspace(0,1000,1000);
%plot(t_axis,light_function(t_axis,chosen_pixel_x,chosen_pixel_y))


for i=1:length(events_x)
   
    if (events_x(i) == chosen_pixel_x)
        if (events_y(i) == chosen_pixel_y)
            chosen_pixel_events_timestamps = [chosen_pixel_events_timestamps events_timestamp(i)];
            chosen_pixel_events_polarity = [chosen_pixel_events_polarity events_polarity(i)];
        end
        
        if (events_y(i) == 0)
            ref_pixel_events_timestamps = [ref_pixel_events_timestamps events_timestamp(i)];
            ref_pixel_events_polarity = [ref_pixel_events_polarity events_polarity(i)];
        end
        
    end 
    
end

%t_axis = linspace(0,1000,1000);


%figure(1)
plot(chosen_pixel_events_timestamps,light_function(chosen_pixel_events_timestamps,chosen_pixel_x,chosen_pixel_y))
figure(2)
plot(ref_pixel_events_timestamps,light_function(ref_pixel_events_timestamps,chosen_pixel_x,0))



%[chosen_freqs,chosen_threshold] = thresh_1_pix(chosen_pixel_x,chosen_pixel_y,chosen_pixel_events_timestamps,chosen_pixel_events_polarity,light_function);


%% Threshold Estimation 

% define light function
light_function = @t_linear;
 
[freqs,threshold] = threshold_estimation(events_x,events_y,events_timestamp,events_polarity,light_function,frame_size_x,frame_size_y);

%% plot frequency and threshold 

y_0 = 125;
x_0 = find(positive_threshold_matrix(y_0,:) == min(positive_threshold_matrix(y_0,:)));
x_1 = find(positive_threshold_matrix(y_0,:) == max(positive_threshold_matrix(y_0,:)));

figure(3)
plot(freqs{y_0,x_0});
hold on
plot(freqs{y_1,x_0});
title('frequency estimation');
legend("c_p = " + positive_threshold_matrix(y_0,x_0) + " c_n = " + negative_threshold_matrix(y_0,x_0),"c_p = " + positive_threshold_matrix(y_0,x_1) + " c_n = " + negative_threshold_matrix(y_0,x_1));
hold off

figure(4)
stem(threshold{y_0,x_0});
title('threshold estimation');
