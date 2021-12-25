% EE5630 HW7 solution: Optimal FIR design, Part I
%
% Yi-Wen Liu
% Dec. 15, 2021
clear; close all;
omega_p = 0.2*pi;
omega_s = 0.24*pi;
x_p = cos(omega_p);
x_s = cos(omega_s);
rng(7414)
Modify = false;         % improve the initializaion way
CompareL = false;       % compare others L
DrawFourMethod = false; % compare others method

Ls = 32; % result will be Type I FIR with 2N+1 samples

% % ---------------------------------------------
% % compare different length
% Ls = [6 12 18];
% CompareL = true; 
% % ---------------------------------------------

% % ---------------------------------------------
% % remove following comment to show comparison
% Ls = 50;
% omega_p = 0.20*pi;
% omega_s = 0.28*pi;
% x_p = cos(omega_p);
% x_s = cos(omega_s);
% DrawFourMethod = true;
% Modify = true;
% % ---------------------------------------------


for ii = 1:length(Ls)

    
    L = Ls(ii);

    figure(2+ii);
    line([-1 x_s],[0 0],'linewidth',1.5); hold on;
    line([x_p, 1],[1 1],'linewidth',1.5);
    set(gca,'ylim',[-0.2 1.2]);
    xlabel('x');

    % Manually select L+2 locations on the x-axis in the ascending order, 
    % with the following constraints:
    % x(1) = -1 and/or x(L+2) = 1
    % one of them, say x(6) = x_s;
    % the next one, say x(7) = x_p;
    df = pi/(L+1);
    ff = pi:-df:0;
    

    xx = cos(ff); xx = xx(:);
    % Change the one closest to x_p to exactly x_p
    % ----------------------original code --------------------------
    [val_p,ind] = min(abs(xx-x_p));
    xx(ind) = x_p;
    if Modify == false
        xx(ind-1) = x_s;
    else
        % --------------------------------------------------------------
        % we further deal with the ind-1-th point to avoid interupting x
        % in transition band
        [val_s,ind_s] = min(abs(xx-x_s)); % find the stopband point
        % it should belong to the right index rather than the ind-1-th
        xx(ind_s) = x_s;        
        % find points in the transition band
        xx_drop = xx(ind_s+1:ind-1);
        % determine whether any point is required to be dropped out.
        while xx_drop ~ []
            % random generate compensated alternating points
            xx_in = 2*rand(length(xx_drop),1) - 1;
            % place it in xx array in sequence
            xx_in = [xx(1:ind_s);xx(ind:end);xx_in];
            xx = sort(xx_in,'ascend');
            % since random generated point may still in the transition band
            % check whether any points are requires to be dropped out once again
            % until suffice the condition
            [val_p,ind] = min(abs(xx-x_p));
            xx(ind) = x_p;
            [val_s,ind_s] = min(abs(xx-x_s));
            xx(ind_s) = x_s;
            xx_drop = xx(ind_s+1:ind-1);
        end
    end

%   xx = WeightingArrange(xx,L,x_p,x_s);
%   xx = NotInTransition(L,omega_p,omega_s);
   
    %% Polynomial fitting
    % Now, find a polynomial p(x) with the following constraints for all x values in the array xx:
    % If x = xx(k) <= x_s then p(x) = (-1)^k * delta, because  it is in the STOPBAND
    % If x = xx(k) >= x_p then p(x) = (-1)^k * delta + 1, because it is in the PASSBAND
    % Reference: Oppenheim & Schafer Eq. (7.113)
    A = zeros(L+2,L+2);
    y = zeros(L+2,1); % y(1:L+1) represents the polynomial, y(end) is delta

    %% STUDENTS: you should set up the A matrix for solving Ay = b
    pow = 0:(L+1);               % power set from 0 to L
    last_column = (-1).^(pow+1); % last column consists of -1 and 1
    A = xx.^pow(1:end-1);        % x_ij^(jj-1)
    A = [A last_column'];        % add last column


    %% setting up the b vector, which represents the desired response
    b = zeros(L+2,1);
    b(ind:end) = 1; % xx(ind) = x_p marks the starting place for the passband.


    y = A\b;

    delta = abs(y(end))
    P = y(1:end-1);

    % Plot the polynomial
    figure(2+ii);
    line([-1 x_s],[delta delta],'linestyle','--','color','red');
    line([x_p 1],[1+delta 1+delta],'linestyle','--','color','red');
    line([-1 x_s],[-delta -delta],'linestyle','--','color','red');
    line([x_p 1],[1-delta 1-delta],'linestyle','--','color','red');


    xplotloc = -1:0.001:1; 
    len = length(xplotloc);
    p_x = zeros(size(xplotloc)); 
    for j = 1:len
        thisterm = 1;
        for m = 1:L+1
            p_x(j) = p_x(j) + P(m)*thisterm;
            thisterm = thisterm * xplotloc(j);
        end
    end
    plot(xplotloc,p_x);
    title(['L = ' num2str(L)])
    %% Find the FIR coefficients
    % P(x) = sum_{k=0:L} a_k x^k = h[0] + 2 sum_{k=1:L} h[k] cos(k arccos(x)) = ...
    %    h[0]+2*sum_{k=1:L} T_k(x), with the summation index going from 1 to L.
    % The next few lines should find (h[0]... h[L]) from 
    % polynomial coeff P(1:L+1)
    h = zeros(2*L+1,1);
    W = myChebyPol(L);  %% STUDENTS: please write this function. 
    % 
    % W = myChebyPol(N) should return an (N+1)x(N+1) matrix, in which the
    % (m+1)th column represents the mth-order Chebyshev polynomial: 
    % T_m(x) = w(1,m+1) + w(2,m+1)x + w(3,m+1)x^2 +... +w(m+1,m+1)x^N for all
    % m = 0 to N.
    % ...
    % ...
    % ...

    %% Converting from the polynomial back to the filter coeffcients
    u = W\P; % the result should be u = [h[0] 2h[1] 2h[2] ... 2h[L]]^T.
    h(1) = u(1); % h[0]
    h(2:L+1) = u(2:end)/2;

    %% Then enforce even symmetry
    h = circshift(h,L);
    h(1:L) = h(end:-1:L+2);

    figure(2+length(Ls)+ii);
    freqz(h);
    title(['L = ' num2str(L)])
    setFontSizeForAll(12)
    [mag, freq] = freqz(h);
    
    if CompareL
        figure(1)
        plot(freq/pi,mag2db(abs(mag)),'DisplayName',['L = ' num2str(L)])
        xlabel('Normalized frequency (\times \pi)')
        ylabel('magnitude response (dB)')
        hold on
        legend
    end
end
if DrawFourMethod
    mag3 = CompareWindowMethod;
    figure
    plot(freq/pi,mag2db(abs(mag3(1,:))),freq/pi,mag2db(abs(mag3(2,:))),freq/pi,mag2db(abs(mag3(3,:))),freq/pi,mag2db(abs(mag)));
    legend('Hann','Hamming','Kaiser','Polynomial')
    xlabel('Normalized frequency (\times \pi)')
    ylabel('magnitude response (dB)')
end
% --------------------------------------------------------------------------
% I even try some other ways, but they are too rubbish.
% function xx = WeightingArrange(xx,L,x_p,x_s)
%     weights = cumsum(cumsum((1:L/2+1)/sum(cumsum(1:L/2+1))));
%     tmp = weights - 1;
%     xx(1:L/2+1) = [-1 tmp(1:end-1)];
%     tmp = 1 - fliplr(weights);
%     xx(L/2+2:end) = [tmp(1:end-1) 1];
%     [~,ind] = min(abs(xx-x_p));
%     xx(ind) = x_p;
%     xx(ind-1) = x_s;
% end
% 
% function xx = NotInTransition(L,omega_p,omega_s)
%     Total_L = pi - omega_s + omega_p;
%     ratio_p = omega_p / Total_L;
%     len_p = ceil(L/2 + 1)
%     len_s = L + 2 - len_p
%     ff1 = linspace(pi,omega_s,len_p);
%     ff2 = linspace(omega_p,0,len_s);
%     ff = [ff1 ff2];
%     xx = cos(ff); xx = xx(:);
% end
% --------------------------------------------------------------------------

function mag = CompareWindowMethod
% 'Hann','Hamming','Kaiser'

for ii = 1:3
    
    sw.winType = ii;
    h = realfilter(sw.winType);
    [mag(ii,:), ~] = freqz(h); % Plot the frequency response
end

function h = realfilter(winType)
    omega_p = 0.2*pi; % edge frequency of passband
    omega_stop = 0.28*pi; % edge frequency of stopband
    A = 65; % In case 3: Kaiser window
    domega = omega_stop - omega_p;
    %% Q1 & Q3
    switch winType
        case 1 % Rectangular window
            N = ceil(2*pi/domega);
        case 2 % Hann window
            N = ceil(4*pi/domega);
        case 3 % Hamming window
            N = ceil(4*pi/domega);
        case 4 % Kaiser window
            % STUDENT: Please estimate beta and the order of Kaiser window 
            % according on A and domega
            if A > 50
                beta = 0.1102*(A-8.7);
            elseif A <= 50 & A >= 21
                beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
            elseif A < 21 & A >= 0
                beta = 0;
            else
                error('A must be equal to zero or a positive value')
            end
            N = ceil((A-8)/(2.285*domega)/2);
            alpha = N;
    end
    % impulse response of filter
    hrecwin = myFIRrecwin((omega_p + omega_stop)/2, N);
    %% Q2 & Q3
    switch winType
        case 1
            h = hrecwin;
        case 2
            win = hann(2*N+1);
            h = hrecwin.*win;
        case 3
            win = hamming(2*N+1);
            h = hrecwin.*win;
        case 4
            % STUDENT: please use besseli() for an explicit expression of w[n]
            % n = (-N:1:N)';
            % if the index of n start from negative N, according to the 
            % document in matlab,
            % win = 1/besseli(0,beta) * besseli(0,beta*sqrt(1-(n/N).^2)); 
            % if the index of n start from 0, according to the lecture note
            n = (0:2*N)';
            win = 1/besseli(0,beta) * besseli(0,beta*sqrt(1-((n-alpha)/alpha).^2));
            h = hrecwin.*win;
    end
    

    function h = myFIRrecwin(omega_c, N)

        n = -N:1:N;
        h = sin(omega_c*n)./(n*pi);
        h(N+1) = omega_c/pi;
        h = h';
    end
end
end