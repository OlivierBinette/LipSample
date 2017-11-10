function pretty_hist( data, varargin )
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
%   pretty_hist(data, [a b])
%       Pretty histogram normalized to be a probability
%       density function on [a, b].
%
%   pretty_hist(data, nbins, [a b])
%       Pretty histogram with _nbins_ bins, normalized to be a probability
%       density function on [a, b].
%
%   pretty_hist(..., 'Type', hist_type)
%       hist_type is one of 'pdf', 'cumulative' and 'counts'. The 'pdf'
%       keyword normalizes the histogram for it to be a probability density
%       function; the 'cumulative' keyword yields a cumulative probability
%       distribution; 'counts' provide a count histogram.
%
%   Examples
%   --------
%   % Basic histogram
%   pretty_hist(rand(1,1000))
%
%   % Use a large number of bins to obtain the empirical cdf.
%   pretty_hist(data, 10000, 'Type', 'cumulative')
%
%   % Counts histogram on [0, 1].
%   pretty_hist(data, [0 1], 'Type', 'counts')
%
%   % Counts histogram on [0, 1] with 50 bins.
%   pretty_hist(data, 50, [0 1], 'Type', 'counts')
%
%   Notes
%   -----
%     - Data points outside of specified limits _[a, b]_ are counted in the
%       outermost bins.
%
%   CC-BY O.B. sept. 15 2017

    % Parsing input arguments.
    
    nstaticin = nargin; % Number of unnamed arguments.
    for i = 1:(nargin-1)
        if ischar(varargin{i})
            nstaticin = i;
            break
        end
    end

    if nstaticin == 1
        a = min(data);
        b = max(data);
        nbins = ceil(3 * numel(data)^0.33);
    elseif nstaticin == 2
        nbins = varargin{1};
        if isequal(size(nbins), [1 1])
            a = min(data);
            b = max(data);
        else
            a = nbins(1);
            b = nbins(2);
            nbins = ceil(3 * numel(data)^0.33);
        end
    elseif nstaticin == 3
        nbins = varargin{1};
        lim = varargin{2};
        a = lim(1);
        b = lim(2);
    end
    
    p = inputParser;
    addOptional(p, 'Type', 'pdf');
    parse(p, varargin{nstaticin:end});
    hist_type = p.Results.Type;
        
    % Drawing the histogram.
    if strcmp(hist_type, 'pdf')
        edges = linspace(a, b, nbins+1) - (b-a)/(2*nbins);
        [N,X] = hist(data, edges(2:end));
         N = nbins * N ./ (sum(N) * (b-a));
    elseif strcmp(hist_type, 'cumulative')
        edges = linspace(a, b, nbins+1) - (b-a)/(2*nbins);
        [N,X] = hist(data, edges(2:end));
        N = cumsum(N ./ sum(N));
    elseif strcmp(hist_type, 'counts')
        edges = linspace(a, b, nbins+1) - (b-a)/(2*nbins);
        [N,X] = hist(data, edges(2:end));
    end
    
    bar(X, N, 1, 'facecolor', [0.8 0.8 0.8], ...
          'edgecolor', [0.8 0.8 0.8]);
    prettify();
    
end
