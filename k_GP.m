% For one-dimensional inputs, computes the covariance function and its partial derivatives wrt:
%  the horizontal lengthscale parameter, l
%  the noise parameter, sigman
%  the vertical lengthscale parameter, sigmaf
%  the frequency, f
%
% N.B. If your data set contains multiple identical x values, rewrite this code to accept the
% indices of x rather than the values themselves; otherwise, the last if-statement below leads
% eventually to the creation of a singular matrix.
%
% Mark Ebden, 2008

function [covar, d_l, d_sigman, d_sigmaf, d_f] = k_GP (x1, x2, theta)

% a) Initializations
if exist('theta')==0,
    theta = [];
end
if length(theta) < 4,
    f = 0;
else
    f = theta(4);
end
if length(theta) < 3,
    sigmaf = 1;
else
    sigmaf = theta(3);
end
if length(theta) < 2,
    sigman = 0;
else
    sigman = theta(2);
end
if length(theta) < 1,
    l = 1;
else
    l = theta(1);
end

% b) Covariance and gradients
covar = sigmaf^2 * exp(-(x1-x2)^2/(2*l^2));
d_l = covar * (l^-3) * (x1-x2)^2; % Differentiate (2.16) from Rasmussen and Williams (2006)
d_sigmaf = 2*sigmaf * exp(-(x1-x2)^2/(2*l^2));
if f > 0,
    covar = covar + exp(-2*sin(pi*f*(x1-x2))^2);
    d_f = exp(-2*sin(pi*f*(x1-x2))^2) * (-4*sin(pi*f*(x1-x2))) * cos(pi*f*(x1-x2)) * f*pi*(x1-x2);
else
    d_f = 0;
end
if x1==x2,
    covar = covar + sigman^2;
    d_sigman = 2*sigman;
else
    d_sigman = 0;
end
