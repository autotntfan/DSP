clc, clear; close all; clear sound;

%{
    DSP HW#9: Echo Cancellation 

    In this homework, we are implementing LMS algorithm.
    Aside from completing the LMS code, you will also have to choose
    a appropriate filter order and a step size.

    Adaptive filter structure:
    desired signal (d): far-end recorded + near-end recorded
    input signal   (x): far-end reference
    output signal  (y): far-end estimated
    error signal   (e): near-end estimated
%}

% parameters
p = ; % filter order 
alpha = ; % step size

% read files
wavFile_d = "double-talk.wav";
wavFile_x = "far-end-reference.wav";
[d, fs_d] = audioread(wavFile_d);
[x_full , fs_x] = audioread(wavFile_x);

% resample wav data to 16k Hz
d = resample(d, 16000, fs_d);
x_full = resample(x_full, 16000, fs_x);

% automatic sync far-end original and double-talk
[x_full, d] = Sync(x_full, d);

% FIR filter
h_pr=[1 -0.9]; 
d = filter(h_pr, 1, d); % pre-emphasis

% AEC initialization
len = length(x_full);
e = zeros(len, 1);
x = zeros(p, 1);
h = zeros(p, 1); % adaptive filter

% LMS Processing
for i = 1 : len
    x = circshift(x, 1);
    x(1) = x_full(i);
    % YOUR CODE BELOW
    y = ;
    e(i) = ; 
    h = h + ; 
    % YOUR CODE ABOVE
end
e = filter(1, h_pr, e); % de-emphasis

% show the final h_hat
figure(1)
plot((1:p)/16000, h);
xlabel('sec');
ylabel('final estim. h(t)');

% write the audio file
audiowrite("LMS.wav", e/max(abs(e)), 16000);

function [x, d] = Sync(x, d)
    min_len = min(length(x), length(d(:,1)));
    d = d(1:min_len, 1);
    x = x(1:min_len);
    [~, xcorr_idx] = max(xcorr(d, x, 'normalized'));
    num_sync_sample = length(x) - xcorr_idx;
    n_delay = abs(num_sync_sample) - 18; 
    x = [zeros(n_delay,1); x];
    len = min(length(d), length(x));
    x = x(1 : len, 1);
    d = d(1 : len, 1);
end