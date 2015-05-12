% This is Matlab code to produce Figure 6 of 'Gaussian Processes for Timeseries Modelling'
% a.k.a. Figures 1 and 2 at www.robots.ox.ac.uk/~mebden/reports/GPtutorial.pdf
% Mark Ebden, Aug 2008.  Comments to mark.ebden@eng.ox.ac.uk
%
% Dependencies: calcGP.m k_GP.m, GPtutorialFcn.m, hypSample.m

%%%%%%%%%%%%%%%%%%%%%%%%% 1. Initialization %%%%%%%%%%%%%%%%%%%%%%%%%

% a) Parameters
sigma_n = 0.3; % Noise variance
f = 0; % f = 1/T.  Set to zero to turn off periodicity
% b) Inputs
X = [-1.5 -1 -.75 -.4 -.25 0]'; n = size(X,1);
y = .55*[-3 -2 -.6 .4 1 1.6]'; % My example
% c) Vector indicating which variables to fit (l, sigma_n, sigma_f, and f)
vars = [1 0 1 0]'; % 0 means "known/given"
% d) Vector of constants for objective function
% (Log transform on the parameter space prevents negative values)
consts = [log([1 sigma_n 1 f])'; vars; n; X; y];

%%%%%%%%%%%%%%%%%%% 2. Fit the ML parameters %%%%%%%%%%%%%%%%%%%%%%%%%

% a) Basic initializations
options = optimset('GradObj','on');
Nstarts = 5;
% b) Choose several varied starting points
l_samp = hypSample ([0.1 10], Nstarts); % for l
sf_samp = hypSample ([0.3 3], Nstarts); % for sigma_f
inits = log([l_samp sf_samp]);
% c) Find candidates
paramVec = [];
for randomStart = 1:Nstarts
    [params, fval] = fminunc(@(params) GPtutorialFcn(params,consts), inits(randomStart,:), options);
    paramVec = [paramVec; fval inits(randomStart,:) params];
end
paramVec(:,2:end) = exp(paramVec(:,2:end));
paramVec(:,1) = paramVec(:,1)/max(abs(paramVec(:,1)));
% d) Select best candidate
paramVec = sortrows(paramVec)
params = paramVec(1,size(inits,2)+2:end);
l = params(1); sigma_f = params(2);

%%%%%%%%%%%%%%%%%%%%% 3. Plot the regression %%%%%%%%%%%%%%%%%%%%%%%%%

% a) Plot Figure 2 in GPtutorial.pdf
xstar = (1:1e3)/1e3 * 2.1 - 1.8;
close all
[xstar, bestEstimate, bounds, K, kstar] = calcGP (@k_GP, X, y, [l sigma_n sigma_f f], 0, xstar, 1);
axis tight, v = axis; axis([-1.7 v(2:4)])
v = axis; line ([v(1) v(1)], [v(3) v(4)], 'Color', 'k')
xlabel ('x'), ylabel('y')
% b) Calculate y*
xindex = 952; % determined by inspection
newX = xstar(xindex); newBounds = bounds(xindex,:);
newEst = bestEstimate(xindex); newSigma = (newBounds(1)-newEst)/1.96;
disp([newEst newSigma.^2]) % Answers reported in Step 2 of page 3 of GPtutorial.pdf
% c) Update Figure 2 with y*
errorbar (newX,newEst,newSigma,newSigma,'Color','green')
plot (newX,newEst,'b.','MarkerSize',15)
errorbar (X,y,sigma_n*ones(size(y)),sigma_n*ones(size(y)),'Color','red', 'Linestyle','none')
plot (xstar,bestEstimate,'k'), plot (X,y,'k.', 'MarkerSize',15)
% d) Plot Figure 1
figure, hold on
errorbar (X,y,sigma_n*ones(size(y)),sigma_n*ones(size(y)),'Color','red','Linestyle','none') % plot error bars
errorbar (newX,newEst,newSigma,newSigma,'Color','green')
text(newX-.01,newEst-newSigma-0.4,'?')
plot (X,y,'k.', 'MarkerSize',15), plot (newX,newEst,'b.', 'MarkerSize',15)
xlabel ('x'), ylabel('y'), axis(v)
