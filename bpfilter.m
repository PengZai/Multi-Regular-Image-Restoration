function H = bpfilter(type, M, N, D0, D1, n)
%bpfilter Computes frequency domain bpfilter filters.
%   H = bpfilter(TYPE, M, N, D0, n) creates the transfer function of
%   a bpfilter filter, H, of the specified TYPE and size (M-by-N).
%   Valid values for TYPE, D0, D1, and n are: 
%
%   'ideal'    Ideal highpass filter with cutoff frequency D0.  n
%              need not be supplied. D0 must be positive.
%
%   'btw'      Butterworth highpass filter of order n, and cutoff
%              D0.  The default value for n is 1.0. D0 must be
%              positive.
%
%   'gaussian' Gaussian highpass filter with cutoff (standard
%              deviation) D0. n need not be supplied. D0 must be
%              positive.



n = 1; % Default value of n.


% Generate bpfilter filter.
lp0 = lpfilter(type, M, N, D0, n);



lp1 = lpfilter(type, M, N, D1, n);
H = lp1 - lp0;