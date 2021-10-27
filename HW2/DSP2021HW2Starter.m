%% EE5630 Digital Signal Processing
% HW2: fundamental frequency (f0) estimation
% Oct. 2021
% Yi-Wen Liu


clear; close all;

DIR = './';
%FILENAME = 'HW2-soprano.m4a';
FILENAME = 'HW2-bass.m4a';

[y,fs1] = audioread([DIR FILENAME]);
y = y(:,1); % choose only the first channel.

%soundsc(y,fs1);
fs = 16000;
p = 20; % linear prediction, no need to change here (beyond the scope of this homework).

y = resample(y,fs,fs1);
%% Parameters to play with
framelen = 0.064; % second. [INVESTIGATE]
L = framelen * fs;

sw.emphasis = 1; % default = 1

numFrames = floor(length(y)/L);

Nfreqs = 2^nextpow2(2*L-1)/2; % Num points for plotting the inverse filter response
df = fs/2/Nfreqs;
ff = 0:df:fs/2-df;

if sw.emphasis == 1
    y_emph = filter([1 -0.95],1,y); % This boosts up the high frequency
                %[PARAM] -0.95 may be tuned anywhere from 0.9 to 0.99
else
    y_emph = y;
end

%% STFT and viewing
%win = ones(L,1); % Rectangular window.
win = hann(L+1); win = win(1:end-1); 
    % Remark: this somewhat awkward way correctly creats a Hann window that
    % satisfies the constant overlap-add condition (will be useful for HW3)
for kk = 1:numFrames
    ind = (kk-1)*L+1:kk*L;
    ywin = y_emph(ind).*win;
    Y = fft(ywin,2*Nfreqs);
    A = lpc(ywin,p); % A direct way to obtain the linear prediction coefficients. 
        
    figure(1);
    subplot(211);
    plot(ind/fs*1000, y(ind));
    xlabel('ms')
    set(gca,'xlim',[kk-1 kk]*framelen*1000);
   
    subplot(212);
    [H,W] = freqz(1,A,Nfreqs);
    Hmag = 20*log10(abs(H));
    Ymag = 20*log10(abs(Y(1:Nfreqs))); %
    Hmax = max(Hmag);
    offset = max(Hmag) - max(Ymag);
    plot(ff,Hmag); hold on;
    plot(ff,Ymag+offset,'r'); hold off;
    set(gca,'xlim',[0 fs/2],'ylim',[Hmax-50, Hmax+5]);
    xlabel('Hz')
    
    drawnow;
    %pause;
end

figure(2)
param.fs = fs;
mySpecgram(y,win,L/2,L,param);

setFontSizeForAll(14); % This function is a plotting routine I created. 
                        % It goes through all the figures and set the font
                        % the same size.