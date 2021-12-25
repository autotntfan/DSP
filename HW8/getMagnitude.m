function mag = getMagnitude(p_x, xplotloc, x)

%% readme
%
% this function helps get the specific magnitude 
% of the polynomial given the x-axis value
%
% INSTRUCTION:
%
% input:    p_x:      polynoimal
%           xplotloc: corresponding x-axis of p_x
%           x:        specific x-axis value(s)
%                     (can be either scalar or vector)
% output:   mag:      magnitude

if max(size(x)) == 1
    mag = p_x(xplotloc >= x);
    mag = mag(1);
else
    mag = zeros(1, max(size(x)));
    for i = 1:max(size(x))
        tmp = p_x(xplotloc >= x(i));
        mag(i) = tmp(1);
    end
end