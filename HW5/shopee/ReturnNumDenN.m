function [NUM, DEN, N] = ReturnNumDenN(omega_c, omega_s, passband, stopband)
    % 根據題目要求的兩個頻率以及兩個頻率響應求出對應的N
    % 根據公式求出N，但因為N要整數因此取了ceil，而這個動作會導致Omega_c*Td或omega變了
    % 我們為了要準確的cutoff freq. 選擇犧牲一些stop freq.的準確度，因為cutoff
    % freq. 仍設定為0.2pi，因此由第一題的式子Td_Omega = 2*tan(omega/2);
    % 知Omega_c*Td仍為同樣的值，但在停止頻率的部分會與第三題得到的稍微不同
    N = ceil(log(((1/stopband)^2-1)/((1/passband)^2-1))/(2*log(tan(omega_s/2)/tan(omega_c/2))));
    Omega_c_Td = 2*tan(omega_c/2)/(1/passband^2-1)^(1/2/N);
    % uniformly distributed angles
    angle = pi/2 + (0.5 : 1 : N-0.5).*pi./N;
    s_p = Omega_c_Td*exp(1i*angle); % look up end of text book
    % map s to z
    z_p = (1+s_p/2)./(1-s_p/2); % Bilinear Transform
    DEN = poly(z_p); % Denominator A(z)

    z_zero = -1*ones(1,N);
    NUM = poly(z_zero); % Numerator B(z)

    DCgain = 2^N/sum(DEN); % B(z)/A(z) evaluated at z = 1, i.e. @ omega = 0.
    NUM = NUM/ DCgain;