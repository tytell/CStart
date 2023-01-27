function [xs,ys] = smooth_position(t, x, y, err)
% SMOOTH_POSITION Smooths x and y coordinates simultaneously.
%   [xs,ys] = smooth_position(t, x, y, err)
%
% Uses `spaps` to produce the smoothest curve that fits through the x and y
% coordinates while maintaining a particular mean squared error. This is
% the suggested algorithm from Walker 1998
% (https://doi.org/10.1242/jeb.201.7.981)
% 
%% Inputs
% *t* Time (sec)
% *x*, *y* x and y positions
% *err* Mean squared error between the smoothed curves and the original.
%       For data in pixels, 0.5 pixels is a good starting point.
%
%% Outputs
% *xs*, *ys* Smoothed x and y positions.
%
%% See also
% smooth_midline

k = find(isfinite(x) & isfinite(y));
a = k(1):k(end);

XY = cat(1, x(:)', y(:)');
sp = spaps(t(k), XY(:,k), err^2*range(t(k)));
XYs = fnval(sp,t(a));

xs = NaN(size(x));
ys = NaN(size(y));

xs(a) = XYs(1,:);
ys(a) = XYs(2,:);
