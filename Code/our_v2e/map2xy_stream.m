function [y_pos, x_pos] = map2xy_stream(map)

% a function that takes a map with the number of events in each pixel and
% converts it to a stream of xy events 

y_pos = [];
x_pos = [];

curr_map = map;

while nnz(curr_map) > 0 
   
    % add one event 
    bin_map = curr_map > 0;
    [y_pos_new,x_pos_new] = find(bin_map);
    
    % remove added events from map 
    curr_map = curr_map - bin_map;
    
    % add new events to stream 
    y_pos = [y_pos ; y_pos_new];
    x_pos = [x_pos ; x_pos_new];
    
end

end 