% The objective function for optimization in GPtutorial.m
% Mark Ebden, 2008

function [fcn, grd] = GPtutorialFcn (params,v)

% a) Initializations
eparams = exp(params); evs = exp(v);
whichVars = v(5:8); % indicators
% Variables:
varCount = 1;
if whichVars(1) == 1,
    l1 = eparams(1);
    varCount = varCount + 1;
end
if whichVars(2) == 1,
    sigma_n = eparams(varCount);
    varCount = varCount + 1;
end
if whichVars(3) == 1,
    sigma_f = eparams(varCount);
    varCount = varCount + 1;
end
if whichVars(4) == 1,
    f = eparams(varCount);
    varCount = varCount + 1;
end
% Constants - whatever's left
if exist('l1') == 0,
    l1 = evs(1);
end
if exist('sigma_n') == 0,
    sigma_n = evs(2);
end
if exist('sigma_f') == 0,
    sigma_f = evs(3);
end
if exist('f') == 0,
    f = evs(4);
end
n = v(9);
X = v(10:n+9); y = v(n+10:end); y = y - mean(y);

% b) Variance and its derivative
K = zeros(n); dKdl = zeros(n); dKds = zeros(n); dKdf = zeros(n); dKdw = zeros(n);
for i = 1:n,
    for j = 1:n,
        [K(i,j), dKdl(i,j), dKds(i,j), dKdf(i,j), dKdw(i,j)] = k_GP (X(i),X(j),[l1 sigma_n sigma_f f]);
    end
end

% c) Calculations
L = chol (K,'lower');
alpha = L'\(L\y);
invK = inv(K);
alpha2 = invK*y; % alpha from page 114, not page 19, of Rasmussen and Williams (2006)

% d) Log marginal likelihood and its gradient
logpyX = -y'*alpha/2 - sum(log(diag(L))) - n*log(2*pi)/2;
dlogp_dl = l1*trace((alpha2*alpha2' - invK)*dKdl)/2;
dlogp_ds = sigma_n*trace((alpha2*alpha2' - invK)*dKds)/2;
dlogp_df = sigma_f*trace((alpha2*alpha2' - invK)*dKdf)/2;
dlogp_dw = f*trace((alpha2*alpha2' - invK)*dKdw)/2;
% NB: Since l = exp(p1), dlogp/dp1 = dlogp/dl dl/dp1 = dlogp/dl exp(p1) = dlogp/dl l, where p1 = params(1)

% e) Format output
fcn = -logpyX; grd = [];
if whichVars(1) == 1,
    grd = [grd -dlogp_dl];
end
if whichVars(2) == 1,
    grd = [grd -dlogp_ds];
end
if whichVars(3) == 1,
    grd = [grd -dlogp_df];
end
if whichVars(4) == 1,
    grd = [grd -dlogp_dw];
end
