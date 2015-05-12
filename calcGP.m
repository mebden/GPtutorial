% Calculates, for xstar x-coordinates, the bestEstimate and the 95%
% confidence interval (bounds) for a Gaussian process regression.
% Mark Ebden, July 2008
%
% Inputs:
%  k - a handle to the covariance function
%  X,y - the training data x- and y-coordinates
%  theta - the GP parameters to be passed to the covariance function
%          (e.g. length, scale, periodicity, etc.
%  sigma_n - the noise (optional -- default 0)
%    NOTE: sigma_n should be left as zero if k *already* includes the noise term, as I usually do.
%  xstar - specify the x-coordinates to solve for (optional -- default is
%          to let the algorithm choose them)
%  graphing - produce a plot of results (optional)
%             Default is to plot the results, unless the function outputs were asked for
%
% Outputs:
%  xstar - x-coordinates (default 100 are chosen)
%  K - covariance matrix among training data
%  V - variance at each xstar

function [xstar, bestEstimate, bounds, K, V] = calcGP (k, X, y, theta, sigma_n, xstar, graphing)

if nargin < 7,
    if nargout == 0,
        graphing = 1; 
    else
        graphing = 0;
    end
end
if nargin < 6 || ~any(xstar),
    r2 =(max(X)-min(X));
    N = 1e3; % How many points to use
    z = .5; % How much to extend the time series on either side (e.g. .5 is 50% extension on left and on right)
    xstar = (0:N)/N * (z*2+1)*r2 + min(X)-r2*z;
end
if nargin < 5,
    sigma_n = 0;
end
if nargin < 4,
    theta = [1 0 1 0];
end

% a) Initializations
meany = mean(y); y = y - meany;
n = size(X,1); lx = length(xstar);
K = zeros(n); kstar = zeros(n,1);
for i = 1:n,
    for j = 1:n,
        K(i,j) = k (X(i),X(j),theta);
    end
end
fstarbar = zeros(lx,1); V = zeros(lx,1);

% b) One-off calculations
diags = max(1e3*eps, sigma_n^2); % beef up the diagonal if sigma_n = 0
L = chol (K + diags*eye(n),'lower'); 
alpha = L'\(L\y);
logpyX = -y'*alpha/2 - sum(log(diag(L))) - n*log(2*pi)/2; % Log marginal likelihood

% c) xstar loop
for q = 1:lx,
    for i = 1:n,
        kstar(i) = k (X(i),xstar(q),theta);
    end
    % Mean of prediction
    fstarbar(q) = kstar' * alpha;
    % Variance of prediction
    v = L\kstar;
    ystar_noise = sigma_n^2; % recall f* has no noise itself (Algorithm 2.1 on p 19 of H152)
    V(q) = k (xstar(q),xstar(q),theta) - v'*v + ystar_noise;
end
bounds = [fstarbar+1.96*sqrt(V) fstarbar-1.96*sqrt(V)]+meany;
bestEstimate = fstarbar+meany;

% d) Output
if graphing == 1,
    figure, hold on
    stdRegion (xstar,bounds)
    plot (xstar,bestEstimate,'k')
    plot (X,y+meany,'b+')
end
