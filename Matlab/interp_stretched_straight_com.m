function [comx, comy] = interp_stretched_straight_com(mx, my, comdist)
% INTERP_STRETCHED_STRAIGHT_COM  Interpolate stretched-straight center of
% mass
%    [comx, comy] = interp_stretched_straight_com(mx, my, comdist)
%
%    Interpolates the x, y location of the stretched straight center of
%    mass (COM), given its distance along the body from the snout.
%
%    mx, my - x and y coordinates of the midline, where columns represent 
%       one time and rows represent one point on the body. The first row 
%       is the head and the last is the tail.
%    comdist - Distance of the stretched-straight COM from the snout
%
%    Returns
%      comx, comy - x and y coordinates of the stretched-straight COM
    arguments
        mx double
        my double
        comdist (1,1) double
    end

% first calculate the arc length along the body
s = [zeros(1, size(mx, 2)); 
    cumsum(sqrt(diff(mx).^2 + diff(my).^2))];

comx = NaN(1, size(mx,2));
comy = NaN(1, size(mx,2));

for fr = 1:size(mx, 2)
    if all(isfinite(mx(:, fr)))
        comx(fr) = interp1(s(:,fr), mx(:,fr), comdist);
        comy(fr) = interp1(s(:,fr), my(:,fr), comdist);
    end
end


        
