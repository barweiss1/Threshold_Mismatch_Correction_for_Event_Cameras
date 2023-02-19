function [events_x,events_y,events_timestamp,events_polarity] = our_v2e(video,gs_frames_array,frames_no)
    

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


end