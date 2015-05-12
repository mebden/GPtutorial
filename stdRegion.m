% Fill a graph with standard deviations; t is nx1 and Range is nx2
% sC is the colour to use, in RGB vector format
% reaxisVar = 0 to leave the axes to be overdrawn, or 1 to redraw them
% Mark Ebden, 2008

function theRange = stdRegion (t, theRange, sC, reaxisVar)
if exist('sC') == 0,
    sC = [1 .8 .8];
end
if exist('reaxisVar') == 0,
    reaxisVar = 0;
end
t = t(:);
fill ([t; flipud(t)], [theRange(:,1); flipud(theRange(:,2))], sC, 'EdgeColor', sC);

if reaxisVar == 1,
    v = axis;
    line (v([1 1]), v([3 4]), 'Color', 'k')
    line (v([1 2]), v([3 3]), 'Color', 'k')
end