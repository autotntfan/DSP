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

clear; close all;

sw.plotFFT = 0;
sw.window = 0; % (4) set 1 for Hann, 0 for Rect

DIR = './';
[x, fs] = audioread([DIR 'boxing.wav']); 
    % ** STUDENTS: find your favorite song and convert it to a .wav file
x = x(1:10*fs); % take the first 10 seconds
x = x(:);
% sound(x,fs)
SigLength = length(x); % signal length
M = 256;               % (3)-f block size

%% The code below designs a low-pass filter using fir1().
%
% type 'help fir1' for more details.
f_cut = 4000; % This Low Pass Filter only reserves parts below 4000Hz
h = fir1(128, f_cut/(fs/2)); % since 0 < Wn < 1, students may use help fir1 to get te idea
h = h(:);
L = length(h); % filter size
%% (1) time cost of conv 
tic
y_conv = conv(x,h);
toc
%% (2) our convolution code
% tic
% L_zp = SigLength + L - 1; % since the final length is L + SigLengh -1, in advance, we append zero
% x_zp = zeros(L_zp,1);
% x_zp(1:SigLength) = x;
% h_zp = zeros(L_zp,1);
% h_zp(1:L) = h;
% y_conv = zeros(L_zp,1);
% inv_x = flip(x_zp);    % according to formulate definition 
% for ii = 1:L_zp
%     y_conv(ii) = sum(h_zp(1:ii).*inv_x(end-(ii-1):end));
% end
% toc
% 
% sum(abs(y_conv - y)) 

%% (3)
% (a)
PDLength = L + M - 1;           % length after convolution
TotalLength = SigLength + L -1; % append zeros such that array size is a power of 2
N_zp = 2^nextpow2(PDLength);    % FFT number of bins to be sufficiently long
% Some preparations before entering blockwise processing
% set length of output
y = zeros(TotalLength, 1);
t = 0:1/fs:TotalLength/fs-1/fs;

% Reserve memory space
h_zp = zeros(N_zp, 1);
h_zp(1:L) = h; % zero-padded impulse response
H = fft(h_zp); % perform FFT on h_zp

if sw.window % Hann window
    win = hann(M+1);
    win = win(1:end-1);
    hopsize = M/2; % determine a hopsize to satisfy the COLA condition 
else % Rectangular window
    win = ones(M,1);
    hopsize = M;
end

num_frame = floor((TotalLength - N_zp)/hopsize) + 1; % total number of frames

% (e) time cost
if ~sw.plotFFT
    tic % start the stopwatch
else
    df = fs/N_zp;
    ff = 0:df:(N_zp-1)*df;   
end

%% Blockwise processing starts here (FFT)
for mm = 1:num_frame
    t_start = (mm-1)*hopsize;
    
    tt = (t_start + 1):(t_start + M);
    x_win = x(tt).*win; % windowing
    x_win_zp = zeros(N_zp, 1);
    x_win_zp(1:M) = x_win;  % zero-padding
    
    % (b) Filtering using direct multiplication in freq. domain
    X =  fft(x_win_zp);       % perform FFT on x_win_zp
    Y =  X.*H;                % multiply the FFT of x_win_zp and h_zp in freq domain
    % (c)
    y_win =  ifft(Y);         % take inverse transform of result in 3.(b)
    % check whether samples beyond 384th are zero
    if sum(y_win(M+L:end)) > 1e-10
        error('wrong signal')
    end
    % (d) overlap-add
    tt2 = 1:N_zp;
    y(t_start + tt2) = y(t_start + tt2) + y_win; % add instead of replacement
    % plot FFT
    if sw.plotFFT
        figure(1);
        Xmag = 20*log10(abs(X));
        Ymag = 20*log10(abs(Y));
        plot(ff/1000, Xmag, 'b', ff/1000, Ymag, 'r--');
        set(gca, 'xlim', [0 (fs/2)/1000], 'ylim', [-100 20]);
        
        xlabel('kHz');
        ylabel('dB');
        title(sprintf('t = %.3f s\n', t_start/fs));
        drawnow;
    end
end

if ~sw.plotFFT
    toc; % end the stopwatch
end


% verify our result
err = abs(y - y_conv);
figure
plot(t,err)
title('error curve')
xlabel('time (s)')
ylabel('error')
sound(y, fs);
