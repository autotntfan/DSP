function [arg, mag] = newArgument(p_x, xplotloc, L, delta, x_p, x_s)
% p_x      : polynomial values at x, p(x) = p0 + p1x + p2x^2 +...+ pLx^L
% xplotloc : x axis
% L        : related to filter length, which is 2L+1
% delta    : tolerance
% x_p      : related to passband frequency and derived from omega_p
% x_s      : related to stopband frequency and derived from omega_s

    % since slopes of two sides near extrema point have different
    % directions
    df = p_x(1:end-1) - p_x(2:end);
    cond = (df(1:end-1) .* df(2:end)) < 0;
    % obtain the indices of extrema points
    Ind_extrema = find(cond==true) + 1;
    % obtain the indices of passband and stopband freq
    [~,Ind_xp] = min(abs(xplotloc-x_p));
    [~,Ind_xs] = min(abs(xplotloc-x_s));
    % arrange indices
    Ind = sort([1 Ind_extrema Ind_xp Ind_xs length(xplotloc)]);
    if length(Ind) < (L+2)
        PlotFigure(xplotloc(Ind),p_x(Ind), length(Ind))
        msg = ['number of alternating points is insufficient, ' ...
               'it requires at least ' num2str(L+2) ' points but gets ' ...
               num2str(length(Ind))];
        error(msg)
    end
    % check where extrema points occur continuously
    while true
        Drop_Ind = [];
        for ii = 2:(length(Ind)-1)
            % deal with alternating points in stopband
            if Ind(ii) <= Ind_xs
                % find the relatively larger one 
                if abs(p_x(Ind(ii-1))) > abs(p_x(Ind(ii)))
                    tmp = ii;
                else
                    tmp = ii - 1;
                end
                % check whether or not adjacent points are invalid
                % if the two points are invalid, then drop the smaller one
                if p_x(Ind(ii-1))*p_x(Ind(ii)) > 0 
                    Drop_Ind = [Drop_Ind tmp];
                end
            % deal with alternating points in passband
            elseif Ind(ii) >= Ind_xp
                if abs(p_x(Ind(ii+1))-1) > abs(p_x(Ind(ii))-1)
                    tmp = ii;
                else
                    tmp = ii + 1;
                end
                if (p_x(Ind(ii+1))-1)*(p_x(Ind(ii))-1) > 0
                    Drop_Ind = [Drop_Ind tmp];
                end
            % alternating points in transitionband is invalid
            else
                Drop_Ind = [Drop_Ind ii];
            end
        end
        % obtain indices
        Drop_Ind = unique(Drop_Ind);
        % drop them
        if ~isempty(Drop_Ind)
            Ind(Drop_Ind) = [];
        else
            break
        end
    end
    % assert length at least equaling to L+2
    if length(Ind) < (L+2)
        PlotFigure(xplotloc(Ind), p_x(Ind), length(Ind))
        msg = ['number of alternating points is insufficient, ' ...
               'it requires at least ' num2str(L+2) ' points but gets ' ...
               num2str(length(Ind))];
        error(msg)
    end
    % check where the start and end point follow rule of delta
    passband_end = (p_x(Ind(end))-1)*(p_x(Ind(end-1))-1);
    stopband_end = p_x(Ind(1))*p_x(Ind(2));
    if passband_end > 0
        Ind = Ind(1:end-1);
    end
    if stopband_end > 0
        Ind = Ind(2:end);
    end
    % if number of points is larger than L+2,then we drop either end or start
    % one which is relatively smaller.
    if length(Ind) > (L+2)
        if abs(p_x(end)-1)>abs(p_x(1))
            Ind = Ind(2:end);
        else

            Ind = Ind(1:end-1);
        end
    end
    % assert length equaling to L+2
    if length(Ind) ~= (L+2)
        PlotFigure(xplotloc(Ind), p_x(Ind), length(Ind))
        error('number of alternating points is wrong')
    end
    arg = xplotloc(Ind);
    mag = p_x(Ind);
    PlotFigure(xplotloc(Ind), p_x(Ind), length(Ind))
    function PlotFigure(arg, mag, len)
        
        figure
        line([-1 x_s],[0 0],'linewidth',1.5); hold on;
        line([x_p, 1],[1 1],'linewidth',1.5);
        line([-1 x_s],[delta delta],'linestyle','--','color','red');
        line([x_p 1],[1+delta 1+delta],'linestyle','--','color','red');
        line([-1 x_s],[-delta -delta],'linestyle','--','color','red');
        line([x_p 1],[1-delta 1-delta],'linestyle','--','color','red');
        set(gca,'ylim',[-0.2 1.2]);
        xlabel('x');
        hold on
        plot(xplotloc,p_x,'b')
        hold on
        plot(arg,mag,'ro')
        txt = join(["There are " num2str(len) " points"]);
        title(txt)
    end
end




    

        
    
    
    
    
    