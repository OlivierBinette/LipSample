function [sample, x, y] = lipsample(f, L, limits, m, varargin)
% Random variates from a Lipschitz continuous probability density function on [a,b].
%
%   s = lipsample(@f, L, [a b], m)
%       Draws _m_ random variates from the probability density _f_ on [_a_, _b_] 
%       which is Lipchitz continuous of order _L_. If _f_ is continuously 
%       differentiable, then the best choice of _L_ is the maximum value 
%       of its derivative.
%
%   s = lipsample(..., 'N', n)
%       ... Uses _n_ mixtures components in the spline envelope of _f_. 
%       The default choice is n = ceil(2*_L_), although increasing _n_ may
%       improve performance in some cases.
%
%       
%   [s, x, y] = lipsample(@f, L, [a b], m)
%       ... Returns the spline envelope constructed by the algorithm: the
%       envelope linearly interpolates the points (x,y).
%
%   Dependencies
%   ------------
%     - Function discretesample.m
%
%   Examples
%   --------
%   % In file myfunc.m
%       function y = myfunc(x)
%           y = 1 + cos(2*pi*x)
%       end
%
%   % A few exact samples
%       sample = lipsample(@myfunc, 2*pi, [0 1], 10000);
%
%   % Plot 10 million variates.
%       sample = lipsample(@myfunc, 2*pi, [0 1], 10000000);
%       hold on
%       pretty_hist(sample, [0 1]);
%       plot(linspace(0,1), myfunc(linspace(0,1)));
%       hold off
%
%   % Plot the envelope constructed by the algorithm
%       [sample, x, y] = lipsample(@myfunc, 4*pi, [0 1], 10000);
%       u = linspace(0, 1, 200);
%       hold on
%       pretty_hist(sample, [0 1]);
%       plot(u, myfunc(u));
%       plot(u, interp1(x,y,u));
%       hold off
%       
%
%   Implementation details
%   ----------------------
%     - Acceptance-rejection sampling. A first degree spline envelope of _f_
%       is constructed. The number of components is a function of _L_, chosen
%       as to maximize expected efficiency.
%
%   Warnings
%   --------
%       _L_ must be greater or equal to the best Lipschitz continuity constant
%       of _f_. Otherwise the algorithm may fail to yield exact samples.
%
%     - Efficiency bottleneck is the evaluation of _f_ at O(m) points. 
%
%   CC-BY O.B. sept. 15 2017

    % Parse input arguments.
    a = limits(1);
    b = limits(2);

    p = inputParser;
    addOptional(p, 'N', ceil(200*L) + 200);
    
    parse(p, varargin{:});
    n = p.Results.N;
        
    % Construct the spline envelope.
    s = (b-a) * L / (2*n);
    x = linspace(0,1,n+1);
    y = arrayfun(f, x*(b-a) + a);
    ylow = arrayfun(f, x*(b-a) + a);
    
    % Use the Lipschitz constant to locally adjust the spline.
    alpha = atan(L);
    d = diff(y);
    beta = abs(atan(n*d/(b-a)));
    r = 0.5*sqrt(((b-a)/n )^2 + d.^2).*sin(pi-alpha-beta)./sin(alpha);
    h = r.*(L - abs(n*d/(b-a)));
    y(1) = y(1) + h(1); ylow(1) = ylow(1) - h(1);
    y(n+1) = y(n+1) + h(n); ylow(n+1) = ylow(n+1) - h(n);
    for i = 2:n
        y(i) = y(i) + max(h(i-1), h(i));
        ylow(i) = ylow(i) - max(h(i-1), h(i));
    end
            
    % Generate random variates following the envelope.
    nProp = ceil((1+s)*m);
    U1 = rand(1, nProp);
    U2 = rand(1, nProp);

    y(1) = y(1)/2;
    y(end) = y(end)/2;
    I = discretesample(y, nProp);
    y(1) = 2*y(1);
    y(end) = 2*y(end);

    U = abs((U1 + U2 + I - 2)/n);
    U(U > 1) = 2 - U(U > 1); % The sample.

    % Generate from  f
    V = rand(1, nProp);
    B = interp1(x, ylow, U);
    passlow = lt(V .* interp1(x,y,U), B);
    sample1 = U(passlow);
    U = U(~passlow); V = V(~passlow);
    sample2 = U(lt(V.*interp1(x,y,U), arrayfun(f, U*(b-a)+a)));
    sample = (b-a)*cat(2, sample1, sample2) + a;
    
    if numel(sample) < m
        sample = cat(2, sample, lipsample(f, L, [a b], m - numel(sample)));
    else
        sample = sample(1:m);
    end
    
    x = x *(b-a) + a;
end