function pretty_hist( data, nbins, a, b)
% Pretty histograms.
%
%   Usage
%   -----
%   pretty_hist(data)
%       Pretty histogram normalized to be a probability density function
%       on [min(data), max(data)].
%
%   pretty_hist(data, nbins)
%       Pretty histogram with _nbins_ bins, normalized to be a probability
%       density function on [min(data), max(data)].
%
%   pretty_hist(data, a, b)
%       Pretty histogram normalized to be a probability
%       density function on [a, b].
%
%   pretty_hist(data, nbins, a, b)
%       Pretty histogram with _nbins_ bins, normalized to be a probability
%       density function on [a, b].
%
%   Examples
%   --------
%   pretty_hist(rand(1,1000))
%
%   Notes
%   -----
%     - Data points outside of specified limits _a_, _b_ are counted in the
%       outermost bins.
%
%   O.B. sept. 15 2017

    % Parsing input arguments.
    if nargin == 1
        a = min(data);
        b = max(data);
        nbins = ceil(3 * numel(data)^0.33);
    elseif nargin == 2
        if isequal(size(nbins), [1 1])
            a = min(data);
            b = max(data);
        else
            a = nbins(1);
            b = nbins(2);
            nbins = ceil(3 * numel(data)^0.33);
        end
    end
        
    % Drawing the histogram.
    edges = linspace(a, b, nbins+1) - (b-a)/(2*nbins);
    [N,X] = hist(data, edges(2:end));
    N = nbins * N ./ (sum(N) * (b-a));
    bar(X, N, 1, 'facecolor', [0.8 0.8 0.8], ...
          'edgecolor', [0.8 0.8 0.8]);
    prettify();
end


