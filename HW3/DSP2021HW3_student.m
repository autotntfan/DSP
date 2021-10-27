% EE5630 DSP
%
% Homework #3 starter code
%
% FIR filter design by windowing, and sequential processing using the 
% Overlap-Add (OLA) analysis and synthesis framework
% 
% Created Oct. 19, 2010. Last updated Oct. 12, 2021.
%
% When turning in the homework, submit your MATLAB code, and a piece of
% music that you use for testing the code.  The music must have human
% vocal parts in it.  

clear; close all; clc;
sw.plotFFT = 0;
sw.window = 1; % set 1 for Hann, 0 for Rect

DIR = './';
[x, fs] = audioread([DIR 'boxing.wav']); 
    % ** STUDENTS: find your favorite song and convert it to a .wav file
x = x(1:10*fs); % take the first 10 seconds
x = x(:);
Ltotal = length(x); % 
M = 128; % 

if sw.window % Hann window
    win = hann(M+1);
    win = win(1:end-1);
    hopsize = ...; % determine a hopsize to satisfy the COLA condition 
else % Rectangular window
    win = ones(M,1);
    hopsize = M;
end

num_frame = ... % total number of frames


%% The code below designs a low-pass filter using fir1().
%
% type 'help fir1' for more details.
f_cut = 4000; % This Low Pass Filter only reserves parts below 4000Hz
h = fir1(128, f_cut/(fs/2)); % since 0 < Wn < 1, students may use help fir1 to get te idea
h = h(:);

%% Some preparations before entering blockwise processing
% set length of output
L = length(h);
y = zeros(length(x) + L - 1, 1);

N_zp = ...;  % Set the FFT number of bins to be sufficiently long

% Reserve memory space
x_win_zp = zeros(N_zp, 1);
h_zp = zeros(N_zp, 1);
h_zp(1:L) = h; % zero-padded impulse response
H = ...; % perform FFT on h_zp

if ~sw.plotFFT
    tic % start the stopwatch
else
    df = fs/N_zp;
    ff = 0:df:(N_zp-1)*df;   
end

%% Blockwise processing starts here
for mm = 1:num_frame
    t_start = (mm-1)*hopsize;
    
    tt = (t_start + 1):(t_start + M);
    x_win = %....; % windowing
    x_win_zp(1:M) = x_win;  % zero-padding
    
    % Filtering using direct multiplication in freq. domain
    X =  ...;       % perform FFT on x_win_zp
    Y =  ...;       % multiply the FFT of x_win_zp and h_zp in freq domain
    y_win =  ...;   % take inverse transform of result in 3.(b)
    
    
    %% overlap-add
    tt2 = 1:(M+L-1);
    y(t_start + tt2) = ...;

    %% plot FFT
    if sw.plotFFT
        figure(1);
        Xmag = 20*log10(abs(X));
        Ymag = 20*log10(abs(Y));
        plot(ff/1000, Xmag, 'b', ff/1000, Ymag, 'r--');
        set(gca, 'xlim', [0 (fs/2)/1000], 'ylim', [-100 20]);
        
        xlabel('kHz');
        ylabel('dB');
        title(sprintf('t = %.3f s\n', t_start/fs));
        %drawnow;
    end


end
if ~sw.plotFFT
    toc; % end the stopwatch
end
sound(y, fs);
