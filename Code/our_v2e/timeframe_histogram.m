function [hist,heatmap] = timeframe_histogram(x,y,ts,pol,Tstart,Tend,Xstart,Xend,Ystart,Yend,nbins,Fignum)

% This function creates a histogram of the timestamps of events in a given
% time frame and space frame
%
% Inputs - 
% x = horizontal pixel index vector of occuring events, i.e x(i) is the x
% index of the i-th event
% y = vetical pixel index vector of occuring events, i.e y(i) is the y
% index of the i-th event
% ts = timestamp vector of occuring events , i.e ts(i) is the timestamp of
% the i-th event
% pol = polarity vector of events 1 is ON event, 0 is OFF event
% Tstart,Tend = start and end of the timeframe 
% Xstart,Xend = start and end of the space frame analyized in x-axis
% Ystart,Yend = start and end of the space frame analyized in y-axis
% Fignum = figure number to plot histogram in
% 
%

% find ts indecies in the relevant time frame
ts_ind = find(ts >= Tstart & ts <= Tend);
% find x,y indecies in the relevant frame
x_ind = find(x >= Xstart & x <= Xend);
y_ind = find(y >= Ystart & y <= Yend);

% merge timeframe conditions with space frame conditions
event_ind = intersect(ts_ind,x_ind);
event_ind = intersect(event_ind,y_ind);

% plot histogram
figure(Fignum);
hist = histogram(ts(event_ind),nbins);
title("Timestamp Histogram between " + Tstart + "[usec] and " + Tend + "[usec]");
xlabel("timestamp[usec]");
ylabel('Number of Events');

X=Xend-Xstart+1;
Y=Yend-Ystart+1;
heatmap=zeros(Y,X);
updated_map=heatmap; % a map to remember what pixel has already triggered

for i=1:length(event_ind)
   
    if(updated_map(Y-(y(event_ind(i))-Ystart),x(event_ind(i))+1-Xstart) == 1)
       continue; 
    end
    heatmap(Y-(y(event_ind(i))-Ystart),x(event_ind(i))+1-Xstart)=ts(event_ind(i))-Tstart;
    updated_map(Y-(y(event_ind(i))-Ystart),x(event_ind(i))+1-Xstart) = 1;
    
end
 


end