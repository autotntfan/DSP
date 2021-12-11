% EE5630 HW6 starter code: FIR filter design by the window method
% Dec 6, 2021
% Prof. Yi-Wen Liu
clear; close all; clc;
%% Parameters
global omega_p omega_stop A
omega_p = 0.2*pi; % edge frequency of passband
omega_stop = 0.28*pi; % edge frequency of stopband
A = 40; % In case 3: Kaiser window
filtername = {'rectangular','Hann','Hamming','Kaiser'};
%% Plotting
for ii = 1:4
    sw.winType = ii;
    figure(ii);
    [n, h] = realfilter(sw.winType);
    stem(n,h); % Plot the impulse response
    xlabel('n')
    ylabel('h')
    title(filtername{ii})
    set(gca,'FontSize',14)
    
    figure(5); 
    hold on
    [mag, w] = freqz(h); % Plot the frequency response
    plot(w/pi,mag2db(abs(mag)),'DisplayName',filtername{ii})
    xlabel('Normalized freq (\times \pi)')
    ylabel('magnitude (dB)')
    title('magnitude response')
    legend('FontSize',20)
    set(gca,'FontSize',16)
    if strcmp(filtername(ii),'Kaiser')
        [h_builtin,w_builtin] = TestBuildInKaiser;
        figure
        plot(w/pi,mag2db(abs(mag)),'black')
        hold on
        plot(w_builtin/pi,mag2db(abs(h_builtin)),'r')
        legend('our','built-in','FontSize',20)
        xlabel('Normalized freq (\times \pi)')
        ylabel('magnitude (dB)')
        title('compare our result with built-in')

    end
    figure(ii+5)
    freqz(h)
    title(filtername{ii})
    set(gca,'FontSize',14)    
end

function [n, h] = realfilter(winType)
    global omega_p omega_stop A
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
                beta = 0.1102*(A-8.7)
            elseif A <= 50 & A >= 21
                beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21)
            elseif A < 21 & A >= 0
                beta = 0
            else
                error('A must be equal to zero or a positive value')
            end
            N = ceil((A-8)/(2.285*domega)/2)
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
    n = 0:1:2*N;
end
function [h, w] = TestBuildInKaiser
    global omega_p omega_stop A
    domega = omega_stop - omega_p;
    N = ceil((A-8)/(2.285*domega)/2);
    hrecwin = myFIRrecwin((omega_stop + omega_p)/2, N);
    if A > 50
        beta = 0.1102*(A-8.7);
    elseif A <= 50 & A >= 21
        beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
    elseif A < 21 & A >= 0
        beta = 0;
    else
        error('A must be equal to zero or a positive value')
    end
    win = kaiser(2*N+1,beta);
    h_kaiser = hrecwin.*win;
    [h, w] = freqz(h_kaiser);
    figure
    stem(0:1:2*N,h_kaiser)
end
