%CELL_FINDER
%
%   [ cells ] = cell_finder( image, pixel_size, border_width,...
%                   minimum_size, maximum_size, eccentricity, solidity,...
%                   max_brightness )
%
%   Takes an image of cells and finds all the cells that it can that
%   conform to a given set of limits. The form of the output is an array of
%   ColiCell objects.
%
%   Cells are filtered by:
%
%   Minimum area - given in nm^2
%   Maximum area - given in nm^2
%
%   Eccentricity (ratio of major to minor axis) - must be greater than 1
%
%   Solidity (proportion of area taken up by cell in convex hull) - between
%           0 and 1.
%
%   Max brightness - maximum mean brightness of cell
%
%   
%   Pixel size gives the size of each pixel, in nm/px (normally on the
%   order of 160
%
%   Border width gives the exclusion border around the edge of the image,
%   in nm
%

function [ cells ] = cell_finder( image_loc, pixel_size, border_width,...
                        minimum_size, maximum_size, eccentricity,...
                        solidity, max_brightness )
                    
load('errors.mat');
                    
% Input validation
% Check for correct number of inputs
error( nargchk( 8, 8, nargin ) );
image = imread( image_loc );
% Ensure image is two dimensional
if ~( numel( image ) > 1 && ndims( image ) == 2 )
    error( 'Image must be two dimensional' )
end

imageSize = size( image );

% Convert filters to correct units
minimumSize = minimum_size / ( pixel_size ^ 2 );
maximumSize = maximum_size / ( pixel_size ^ 2 );
borderWidth = ceil( border_width / pixel_size );

% Find cell edges
cellEdges = edge(...
    imfilter( image, fspecial( 'gaussian', 10, 15 ) ),...
    'canny' );

borderMask = zeros( size( image ) );
borderMask( borderWidth + 1:size( image, 1 ) - borderWidth - 1,...
    borderWidth + 1:size( image, 2 ) - borderWidth - 1 ) = 1;
%{
cellsMask = imopen(...
                bwmorph(...
                    imfill(...
                        bwmorph(...
                            bwmorph( cellEdges & borderMask, 'thicken' ),...
                        'close' ),...
                    'holes' ),...
                'thin' ),...
            strel( 'disk', 4 ) );
%}
cellsMask = imopen(...
                imfill(...
                    bwmorph( cellEdges & borderMask, 'close' ),...
                'holes' ),...
            strel( 'disk', 4 ) );

checks = { 'Image', 'BoundingBox', 'MeanIntensity' };
if check_arg( minimum_size ) || check_arg( maximum_size )
    checks( size( checks, 2 ) + 1 ) = { 'Area' };
end

if check_arg( solidity )
    checks( size( checks, 2 ) + 1 ) = { 'Solidity' };
end

if check_arg( eccentricity )
    checks( size( checks, 2 ) + 1 ) = { 'MajorAxisLength' };
    checks( size( checks, 2 ) + 1 ) = { 'MinorAxisLength' };
end

rp = regionprops( cellsMask, image, checks );

cells = cell( 1, numel( rp ) );

for i = 1:numel( rp )
    valid = 0;
    
    if check_arg( minimum_size )
        if rp(i).Area < minimumSize
            valid = TOO_SMALL;
        end
    end
    
    if check_arg( maximum_size )
        if rp(i).Area > maximumSize
            valid = TOO_BIG;
        end
    end
    
    if check_arg( eccentricity )
        if rp(i).MajorAxisLength / rp(i).MinorAxisLength < eccentricity
            valid = TOO_ROUND;
        end
    end
    
    if check_arg( solidity )
        if rp(i).Solidity < solidity
            valid = TOO_BENT;
        end
    end
    
    if check_arg( max_brightness )
        if rp(i).MeanIntensity > max_brightness
            valid = TOO_BRIGHT;
        end
    end
    
    cells{i} = ColiCell( rp(i).Image, rp(i).BoundingBox(1),...
        rp(i).BoundingBox(2), rp(i).MeanIntensity, valid, imageSize );
    
end


end



function ok = check_arg( arg )

ok = isa( arg, 'double' ) && numel( arg ) == 1 && arg > 0;

end