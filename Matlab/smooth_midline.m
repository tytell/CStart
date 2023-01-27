function [mxs,mys, dxds,dyds] = smooth_midline(t,mx,my,len, options)
% SMOOTH_MIDLINE  Smooths all of the points in a fish midline simultaneously
%    function [mxs,mys, dxds,dyds] = smooth_midline(t,mx,my,len, options)
%    
%    Smooths x and y coordinates of a fish midline simultaneously in space
%    (along the body) and in time. Uses `spaps` to find the smoothest
%    surface that passes through all of the body points over time with a
%    particular spatial and temporal mean squared errors.
%
%    t - Time (sec)
%    mx, my - Coordinates of the midline, where columns represent one time
%       and rows represent one point on the body. The first row is the head
%       and the last is the tail.
%    len - Fish body length
%    options.serr - Mean squared error across the body in space. For
%       measurements in pixels, 0.5 pix is a good starting value.
%    options.terr - Mean squared error across each point in time.
%    options.derivatives - Return the spatial derivatives of the midline
%        (true or false).
%
%    Returns
%      mxs,mys - Smoothed x and y coordinates of the midline
%      dxds,dyds - Spatial derivatives of x and y
%    
%    See also smooth_position
arguments
    t (1, :) double
    mx (:, :) double
    my (:, :) double
    len (1,1) double
    options.serr (1,1) double = 0.5
    options.terr (1,1) double = 0.5
    options.derivatives (1,1) logical = false
end

nfr = size(mx,2);
npt = size(mx,1);

% this will only work if all of the points are not NaN
s = [zeros(1,nfr); cumsum(sqrt(diff(mx).^2 + diff(my).^2))];

warning('smoothMidline:lengthinfo', ...
    'Digitized length %.2f+-%.2f = %.1f%% of input (%.2f)\n', ...
    mean(s(end,:), "omitnan"), std(s(end,:), "omitnan")/sqrt(sum(isfinite(s(end,:)))), ...
    mean(s(end,:), "omitnan")/len*100, len);

if ((abs((mean(s(end,:), "omitnan")-len)/len) > 0.05) || ...
        (std(s(end,:), "omitnan")/len  > 0.05))
    warning("smoothMidline:lengthvaries",'Lengths don''t match up (either digitized length varies or is different from input).');
end

serr = options.serr;
terr = options.serr;

ds = diff(s);
ds0 = len/(npt-1);
if (any((ds - ds0)/ds0 > 0.05))
    error("smoothMidline:pointorder",'Distance along midline increases weirdly.');
end

s0 = median(s,2, "omitnan");
s = linspace(0,s0(end),npt);

k = find(all(isfinite(mx) & isfinite(my)));
a = k(1):k(end);

XY = cat(1,shiftdim(mx(:,k),-1),shiftdim(my(:,k),-1));
sp = spaps({s0,t(k)}, XY, ...
    {serr^2*range(s0)*length(t(k)) ...
    terr^2*range(t(k))*length(s0)});
XYs = fnval(sp,{s,t(a)});

mxs = NaN(size(mx));
mys = NaN(size(my));
mxs(:,a) = squeeze(XYs(1,:,:));
mys(:,a) = squeeze(XYs(2,:,:));

if options.derivatives
    dxyds = fnval(fnder(sp,[1 0]),{s,t(a)});
    dxds = squeeze(dxyds(1,:,:));
    dyds = squeeze(dxyds(2,:,:));
else
    dxds = [];
    dyds = [];
end

end