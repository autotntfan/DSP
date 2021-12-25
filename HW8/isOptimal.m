function IsOptimal = isOptimal(mag, delta)

thres = 0.05*delta; % 0.005

len = length(mag);
error = zeros(1,len);
for i = 1:len
    if mag(i) < 0.5 % the extremum is in stopband
        error(i) = abs(mag(i)) - delta;
    elseif mag(i) > 0.5 % the extremum is in passband
        if mag(i) > 1
            error(i) = mag(i) - (1+delta);
        elseif mag(i) < 1
            error(i) = (1-delta) - mag(i);
        end
    end
end
   
check = error > thres;
IsOptimal = all(~check); 