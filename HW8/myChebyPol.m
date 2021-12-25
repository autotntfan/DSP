function  W = myChebyPol(N)
% each column of W represents the n-th Tn(x)
% each row represents the coefficients of x^0,x^1,...,x^L
% e.g. T1(x) = 1, so the first column is [1,0,0,...,0]^T
%      T2(x) = x, so the second column is [0,1,0,...,0]^T
    W = zeros(N+1,N+1);
    % T0(x) term
    T0 = zeros(N+1,1);
    T0(1) = 1;
    % place T0(x) term into the first column of W
    W(:,1) = T0;
    % T1(x) term
    T1 = zeros(N+1,1);
    T1(2) = 1;
    % place T1(x) term into W
    W(:,2) = T1;
    % coef1 is 2x term
    coef1 = zeros(N+1,1);
    coef1(2) = 2;
    % coef2 is -1 term
    % Tn(x) = 2xTn(x-1) - 1Tn(x-2);
    coef2 = zeros(N+1,1);
    coef2(1) = -1;
    for ii = 3:(N+1)
        tmp = conv(W(:,ii-1),coef1) + conv(W(:,ii-2),coef2);
        W(:,ii) = tmp(1:N+1);
    end
end
