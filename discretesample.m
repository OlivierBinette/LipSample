function k = discretesample( p, m )
%  Samples from a discrete distribution on the integers.
%   
%   k = discretesample(p, m)
%       Draws _m_ samples from the discrete probability distribution on
%       (1, 2, ..., numel(p)) specified by _p_.
%
%   Examples
%   --------
%   % 10 samples from the uniform distribution on (1,2,3,4,5).
%   p = ones(5)/5;
%   sample = discretesample(p, 10);
%
%   % Histogram of 10 000 samples.
%   p = 1:100;
%   p = p/sum(p);
%   hist(discretesample(p, 10000);
%
%   Warnings
%   --------
%     - Not compatible with some versions of Octave.
%
%   O.B. sept. 15 2017

    % Parameter (partial) validation.
    s = sum(p);
    if s == 0
        error('Probabilities sum to 0.')
    end

    % Sampling.
    c = cumsum(p);
    edges = [0 c];
    [~, k] = histc(rand(1, m)*c(end), edges);
end