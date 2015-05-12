% Sample a hyperbola N times for a parameter bounded by the two components of 'bounds'.
% Example usage: if p(x) = k/x from 0.1 to 4, and p(x) = 0 otherwise:
%                x = hypSample ([.1 4], 1e5); hist(x,100)
% Bounds can be calculated as per Section 5.5.1 in 'Data Analysis: A Bayesian Tutorial',
% 2nd edition, D.S. Sivia and J. Skilling, 2006.
%
% Mark Ebden, July 2008

function x = hypSample (bounds, N)

xmin = bounds(1); xmax = bounds(2);
F = rand(N,1);
x = xmin.^(1-F) .* xmax.^F;
