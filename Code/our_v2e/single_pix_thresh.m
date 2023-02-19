function [m,rms] = single_pix_thresh(ts,pol,T,plot_flag)

% This function is a supporting function of OffE_thresh_est for the threshold estimation in
% OffE method from "Event Camera Cabliration of Per-pixel Biased Contrast
% Threshold". This function operates on a single pixel and applies curve
% fitting of the sum of polarities to the total event number and calculate
% the RMS of the estimation. 
%
% Inputs - 
% ts = timestamp vector of occuring events for a single pixel
% pol = polarity vector of events 1 is ON event, 0 is OFF event
% T = vector that contains time domain division to estimation segments
% plot_flag = debug flag that plot the estimation graph
% 
% Outputs - 
% m = spline of the linear curve fitting 
% rms = RMS of the estimation to the data 
% 
%

% convert polarities to {-1,1} representation 
signed_pol = (-1).^(pol+1);

% define curve fitting parameters
seg_num = length(T);
event_num = zeros(seg_num-1,1); % sum(abs(sigma_i)) in the paper
pol_sum = zeros(seg_num-1,1); % sum(sigma_i) in the paper

% calculate event number and polarity sum at each window 
for i = 1:length(T)-1
   
    % find events matching to the timing window
    idx = find( ts >= T(i) & ts < T(i+1) );
    event_num(i) = length(idx);
    pol_sum(i) = sum(signed_pol(idx));
    
end

% linear curve fit to the model - pol_sum = m*event_num + d
A = [event_num ones(size(event_num))];

if rcond(A'*A) < 1e-10
    rms = NaN;
    m = NaN;
else

    vec = (A'*A)\A'*pol_sum;
    m = vec(1);

    % calculate RMS
    rms = sqrt(mean((polyval(vec,event_num) - pol_sum).^2));


end
% plot estimation graph for debug 
if plot_flag == true
    
    fit = polyval(vec,event_num);
    figure;
    plot(event_num,fit);
    hold on
    plot(event_num,pol_sum,'o');
    hold off
    title('Polarity Sum Vs. Event Num Curve Fit');
    xlabel('Event Number');
    ylabel('Polarity Sum');
    legend('Curve Fit','Data');
    
end


end