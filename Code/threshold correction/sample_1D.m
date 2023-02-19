function [ts,pol] = sample_1D(x,t,thres_pos,thres_neg)

% this function recieves a signal, a time vector and a threshold value,
% samples it like an event camera - providing polarities and and timestamps at the times the threshold was crossed

%  ----- code -----

delta = x(1); % this will remember the accumulating change in x for sampling 

% initialize ts and pol vectors
ts = [];
pol = [];

len = length(x);
x = [x 0]; % add zero at the end so no sample is lost

% run on x and add differences until the threshold is met
for n = 1:len
    
    % calculate polarity and select correct threshold
    curr_pol = sign(delta);
    pol_01 = 0.5*(curr_pol+1); % convert to {0,1} representation instread of {+1,-1}
    thres = pol_01*thres_pos + (1-pol_01)*thres_neg; 
    
    % calculate how many events where generated at this timestamp
    events_num = floor(abs(delta/thres));
    
    % add event to ts and pol vectors
    while events_num > 0
       ts = [ts t(n)];
       pol = [pol curr_pol];
       events_num = events_num-1;
       delta = delta - thres;
    end
    
    % add accumulated difference for next sample
    delta = delta + x(n+1) - x(n);
    
end



end