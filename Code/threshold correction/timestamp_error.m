function [fst_sigma_pre,fst_sigma_post,lst_sigma_pre,lst_sigma_post] = timestamp_error(ts_pre,ts_post,Tstart,Tend)


% reconstruct the signal simply by incresing or decreasing the amplitude by
% the threshold value every time an event occurs

%  ----- code -----

% get image dimensions 
[pixel_x,pixel_y] = size(ts_pre);

% initialize matrices for first and last timestamp
t_first_pre = zeros(size(ts_pre));
t_last_pre = zeros(size(ts_pre));

t_first_post = zeros(size(ts_pre));
t_last_post = zeros(size(ts_pre));

% run over all pixels 
for x = 1:pixel_x
   for y = 1:pixel_y
       
       % find pre correction first and last timestamps
       t_first_pre(x,y) = ts_pre{x,y}(find( ts_pre{x,y} > Tstart & ts_pre{x,y} < Tend,1,'first'));
       t_last_pre(x,y) =  ts_pre{x,y}(find( ts_pre{x,y} > Tstart & ts_pre{x,y} < Tend,1,'last'));
              
       % find pre correction first and last timestamps
       t_first_post(x,y) =  ts_post{x,y}(find( ts_post{x,y} > Tstart & ts_post{x,y} < Tend,1,'first'));
       t_last_post(x,y) =  ts_post{x,y}(find( ts_post{x,y} > Tstart & ts_post{x,y} < Tend,1,'last'));
       
   end
end

% calculate std of the timestamps for all pixels 
fst_sigma_pre = std(t_first_pre(:));
fst_sigma_post = std(t_first_post(:));
lst_sigma_pre = std(t_last_pre(:));
lst_sigma_post = std(t_last_post(:));


end