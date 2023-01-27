function [comx, comy] = get_com(mx, my, options)
% GET_COM  Estimates the location of the center of mass of a fish
%    [comx, comy] = get_com(mx, my, options)
% 
%    Estimates the location of the center of mass (COM) for a fish, given
%    the coordinates of the midline and the body width and height, or the
%    mass of each segment of the body.
%
%    mx, my - x and y coordinates of the midline, where columns represent one time
%       and rows represent one point on the body. The first row is the head
%       and the last is the tail.
%    options.width - Width (horizontal) of the body, in fractions of the
%       body length (so it should be between 0 and about 0.5)
%    options.height - Height (vertical) of the body in fractions of body
%       length (optional)
%    options.segmentmass - Mass of each segment along the body.
%
%    You must pass in options.width or options.segmentmass. If you also
%    give options.height, it will calculate the mass assuming that the body
%    width and height vary and the body is constant density.
%
%    Returns
%      comx, comy - x and y locations of the center of mass
    arguments
        mx double
        my double

        options.segmentmass = []
        options.width = []
        options.height = []
    end
    npt = size(mx, 1);
    nfr = size(mx, 2);
    
    if isempty(options.segmentmass) && isempty(options.width)
        error('Requires at least the width of the body or the mass per length as a function of length');
    end
    
    if any(options.width > 0.5) || any(options.height > 1)
        warning('width and height should be in fractions of body length.')
    end
    
    if isempty(options.segmentmass)
        s0 = linspace(0, 1, npt)';
        if (length(options.width) ~= npt)
            warning('Width has a different number of points than midline. Interpolating, assuming equal spacing');
            s1 = linspace(0, 1, length(options.width));
        
            options.width = interp1(s1, options.width, s0);
        end
    
        if isempty(options.height)
            options.height = repmat(mean(options.width), [length(options.width) 1]);
        elseif length(options.height) == 1
            options.height = options.height([1 1]);
        end
    
        if (length(options.height) ~= npt)
            warning('Height has a different number of points than midline. Interpolating, assuming equal spacing');
            s1 = linspace(0, 1, length(options.height));
        
            options.height = interp1(s1, options.height, s0);
        end
    
        % width and height should be the total width and height, but the
        % formulas below need the radii
        options.width = options.width / 2;
        options.height = options.height / 2;
        
        dw = diff(options.width);
        dh = diff(options.height);
        dl = diff(s0);
    
        options.width = options.width(1:end-1);
        options.height = options.height(1:end-1);
    
        % we treat each segment as a truncated cone with an oval-cross section
        % and width and height that depend on length
        % this is actually just the segment volume. We're assuming constant
        % density so that it's proportional to the mass
        options.segmentmass = pi * (options.width.*options.height.*dl + ...
            0.5 * (dw./dl .* options.height + dh./dl .* options.width) .* dl.^2 + ...
            1/3 * dw .* dh .* dl);
    elseif length(options.segmentmass) ~= npt - 1
        error('Length of segmentmass should be one less than the number of body points');
    end
    
    % first calculate the arc length along the body
    s0 = [zeros(1, size(mx, 2)); 
        cumsum(sqrt(diff(mx).^2 + diff(my).^2))];
    % midpoint of each segment
    s1 = (s0(1:end-1,:) + s0(2:end,:))/2;
    
    % interpolate the location of each midpoint
    mx1 = NaN(size(s1));
    my1 = NaN(size(s1));
    for fr = 1:nfr
        if all(isfinite(s0(:,fr))) && all(isfinite(mx(:,fr))) && ...
                all(isfinite(my(:,fr)))
            mx1(:,fr) = interp1(s0(:,fr), mx(:,fr), s1(:,fr));
            my1(:,fr) = interp1(s0(:,fr), my(:,fr), s1(:,fr));
        end
    end
    
    totalmass = sum(options.segmentmass);
    segmentmass = repmat(options.segmentmass(:), [1 nfr]);
    
    % and estimate the COM location
    comx = sum(mx1 .* segmentmass) / totalmass;
    comy = sum(my1 .* segmentmass) / totalmass;
end    

