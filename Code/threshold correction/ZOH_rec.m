function x_zoh = ZOH_rec(ts,pol,thres_pos,thres_neg,t)

% reconstruct the signal simply by incresing or decreasing the amplitude by
% the threshold value every time an event occurs

%  ----- code -----

x_zoh = zeros(1,length(t));
Ts = t(2) - t(1);

for i = 1:length(ts)

    if pol(i) == 1
        x_zoh = x_zoh + thres_pos*heaviside((t-ts(i))/Ts);
    else
        x_zoh = x_zoh + thres_neg*heaviside((t-ts(i))/Ts);
    end
   
    
end