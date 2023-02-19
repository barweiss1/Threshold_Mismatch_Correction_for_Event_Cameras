function ll = linlog(y)

% implement lin-log function which is linear for y<=20 and
% logarithmic for y>20
% input image should be represented in uint8 format
func = @linlog_single;
ll = arrayfun(func,y);

end


function l = linlog_single(y)

if y<=20
    l = log(20)*y/20;
elseif y>20
    l = log(y);
end

end