function [sample, x, y] = lipsample(f, L, limits, m, varargin)
% Random variates from a Lipschitz continuous probability density function on [a,b].
%
%   s = lipsample(@f, L, [a b], m)
%       Draws _m_ samples from the probability density _f_ on [_a_, _b_] 
%       which is Lipchitz continuous of order _L_. If _f_ is continuously 
%       differentiable, then the best choice of _L_ is the maximum value 
%       of its derivative.
%
%   s = lipsample(..., 'N', n)
%       ... Uses _n_ mixtures components in the spline envelope of _f_. 
%       The default choice is n = ceil(2*_L_), although increasing _n_ may
%       improve performance in some cases.
%
%   s = lipsample(..., m, 'Tolerance', epsilon)
%       ... When _m_ is large, samples follow an approximation of _f_
%       that is at distance at most epsilon in the supremum norm. By
%       default, epsilon = 0.0005. If epsilon = 0, exact samples are drawn
%       for all values of m.
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
%   % Plot 10 million samples based on a low tolerance approximation of _f_.
%       sample = lipsample(@myfunc, 2*pi, [0 1], 10000000, 'Tolerance', 0.0001);
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
    addOptional(p, 'Tolerance', 0.0005);
    addOptional(p, 'N', ceil(50*sqrt(L)) + 1);
    
    parse(p, varargin{:});
    tolerance = p.Results.Tolerance;
    n = p.Results.N;
        
    % Construct the spline envelope.
    s = (b-a) * L / (2*n);
    x = linspace(0,1,n+1);
    y = arrayfun(f, x*(b-a) + a);
    
    % Use the Lipschitz constant to locally adjust the spline.
    M = abs(atan(n*diff(y)/(b-a)));
    r = sqrt(((b-a)/n)^2 + diff(y).^2);
    alpha = atan(L);
    h = ((b-a)*tan(alpha)/n - r.*sin(M))/2;
    y(1) = y(1) + h(1);
    y(n+1) = y(n+1) + h(n);
    for i = 2:n
        y(i) = y(i) + max(h(i-1), h(i));
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
    
    % Evaluate f at the points U*(b-a)+a
    if (tolerance > 0) && (m > 0.5*(b-a) * L / tolerance) % Optional approximation when m is large
        nDiscretize = ceil(0.5 * (b-a) * L / tolerance)+2;
        fx = linspace(a, b, nDiscretize);
        fy = arrayfun(f, fx);
        B = interp1(fx, fy, U*(b-a) + a);
    else
        B = arrayfun(f,U*(b-a) + a);
    end
    
    % Generate from  f
    V = rand(1, nProp);
    sample = U(lt(V .* interp1(x,y,U), B))*(b-a)+a;
    
    if numel(sample) < m
        sample = cat(2, sample, lipsample(f, L, [a b], m - numel(sample)));
    else
        sample = sample(1:m);
    end
    
    x = x *(b-a) + a;
end