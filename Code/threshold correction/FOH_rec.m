function x_foh = FOH_rec(ts,pol,thres_pos,thres_neg,t)

% reconstruct the signal simply by incresing or decreasing the amplitude by
% the threshold value every time an event occurs

%  ----- code -----

x_foh = zeros(1,length(t));

[ts_uniq,pol_uniq,e_num] = uniqify_ts_pol(ts,pol);

pol_01 = 0.5*(pol_uniq+1); % convert to {0,1} representation instread of {+1,-1}
thres = pol_01*thres_pos + (1-pol_01)*thres_neg; % convert polarities to threshold for each event

% add inital value if the first sample exists
if t(1) == ts(1)
   x_foh(1) = e_num(1)*thres(1);
end

% add first timestamp 
end_idx = find( t == ts_uniq(1) );
idx = 1:end_idx;
if e_num(1) > 1 
    % if there are several events at the same timestamp than use a step
    % rise rather than a linear one since that indicates a rise faster
    % than the temporal sample rate
    linear = e_num(1)*thres(1)*ones(1,length(idx)-1);
    x_foh(idx(2:end)) = linear + x_foh(idx(1));
else
    linear = linspace(0,thres(1),length(idx)-1);
    x_foh(idx(2:end)) = linear + x_foh(idx(1));
end


% add the rest fo the timestamps
for i = 2:length(ts_uniq)
    
    % find relevant indecies
    idx = find( t >= ts_uniq(i-1) & t <= ts_uniq(i));
    
    if e_num(i) > 1 
        % if there are several events at the same timestamp than use a step
        % rise rather than a linear one since that indicates a rise faster
        % than the temporal sample rate
        linear = e_num(i)*thres(i)*ones(1,length(idx)-1);
    % if the change is over a long amout of time the linearization won't
    % perform well so we use stepwise approximation
    elseif ts_uniq(i) - ts_uniq(i-1) >= 0.2
        linear = zeros(1,length(idx)-1);
        linear(end) = thres(i);
    else
        linear = linspace(0,thres(i),length(idx)-1);
    end
    x_foh(idx(2:end))  = linear + x_foh(idx(1));
    
end

x_foh(idx(end):end) = x_foh(idx(end));

end
