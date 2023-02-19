function [ts_uniq,pol_uniq,e_num] = uniqify_ts_pol(ts,pol)

% function that counts same timestamps and counts how many events occured
% at every time 

ts_uniq = [];
pol_uniq = [];
e_num = [];

len = length(ts);
i = 0;

while i < len ; i = i+1;

    % add time stamp to the list 
    ts_uniq = [ts_uniq ts(i)];
    pol_uniq = [pol_uniq pol(i)];
    e_num = [e_num 1];
    
    % count similar events
    while i < len && ts(i) == ts(i+1)
        % add event count
        e_num(end) = e_num(end)+1;
        i = i+1;
       
    end
    
end