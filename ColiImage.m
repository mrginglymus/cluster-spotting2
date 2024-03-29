%COLIIMAGE
%
%   This class handles a pair of images, one containing CheZ information
%   and one containing CheY information. It contains a function to locate
%   cells in the image, and a function to locate clusters in each cell,
%   along with functions to plot all of these data.
%


%   Class must inherit handle so it can be self modifying
classdef ColiImage < handle
    
    properties
        
        cellImage           %   location of the cell image
        clusterImage        %   location of the cluster image
        imageSize           %   size of the image
        cells               %   cell array of each cell
        pixelSize           %   pixel size of image, rather than camera
                            %       i.e., 16�m camera w/ 100X objective =
                            %       160nm pixel size
        
    end
    
    methods
        
        %   Constructor function. Takes cell and cluster image locations
        %   and the pizel size
        function self = ColiImage( cellImageLocation, clusterImageLocation,...
                pixelSize )
            
            %   Calculate and store the important values
            self.cellImage = cellImageLocation;
            self.clusterImage = clusterImageLocation;
            self.imageSize = size( imread( self.cellImage ) );
            self.pixelSize = pixelSize;
            
            self.getCells( 1000, 5000000, 15000000, 1.3, 0.85, [], 500 );
            
            self.getClusters( 1.5, 600, 250, 750 );
                                  
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %   getCells - function to find and categorise every cell in the  %
        %   image.                                                        %  
        %                                                                 %
        %   Inputs:                                                       %
        %       border      -   exclusion zone around edge of image to    %
        %                       prevent picking up half cells             %
        %       minSize     -   minimum size of cell to be considered, in %
        %                       square microns                            %
        %       maxSize     -   as above                                  %
        %       aspectRatio -   minimum ratio of major axis length to     %
        %                       minor axis length of cell                 %
        %       solidity    -   minimum ratio of cell area to convex area %
        %                       of cell; removes bent/joined cells        %
        %       maxBright   -   maximum mean brightness of cell           %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function getCells( self, border, minSize, maxSize, aspectRatio,...
                solidity, maxBright, expansion )
            
            %   load the error definition file
            load( 'errors.mat' );
            
            %   read in the cell image
            image = imread( self.cellImage );
            
            %   convert constraints from microns to pixels
            minSize = minSize / ( self.pixelSize ^ 2 );
            maxSize = maxSize / ( self.pixelSize ^ 2 );
            border = ceil( border / self.pixelSize );
            
            expansion = ceil( expansion / self.pixelSize );
            
            %   get the black level
            blackLevel = self.getBlackLevel( 'cell', 16 );
            
            %   locate the cell edges using the canny filter. This can only
            %   tolerate a moderate level of noise, so run a quick gaussian
            %   over it first
            cellEdges = edge(...
                imfilter( image, fspecial( 'gaussian', 10, 15 ) ),...
                    'canny' );
            
            %   create a border mask to eliminate any cells touching the
            %   edge
            borderMask = zeros( self.imageSize );
            borderMask( border + 1:self.imageSize(1) - border - 1,...
                border + 1:self.imageSize(2) - border - 1 ) = 1;
            
            %   bwmorph - close any small gaps in the edges of the cells
            %   imfill  - fill in the cell outlines
            %   imopen  - remove any small blobs that are clearly not cells
            cellMask = imopen(...
                          imfill(...
                              bwmorph( cellEdges & borderMask, 'close' ),...
                          'holes' ),...
                       strel( 'disk', 4 ) );
            
            %   set up which checks to make on the cells
            checks = { 'Image', 'BoundingBox', 'MeanIntensity', 'Area',...
                'Solidity', 'MajorAxisLength', 'MinorAxisLength' };
            
            %   calculate the metrics of each cell
            rp = regionprops( cellMask, image, checks );
            
            %   calculate mean intensity of CheZ
            rpz = regionprops( cellMask, imread( self.clusterImage ),...
                'MeanIntensity' );
            
            %   prepare the cellarray for filing
            self.cells = cell( 1, numel( rp ) );
            
            %   loop through each cell
            for i = 1:numel( rp )
                %   monitor the status of the cell
                valid = 0;
                
                %   check the cell is not too small, if requested
                if check_arg( minSize )
                    if rp(i).Area < minSize
                        valid = TOO_SMALL;
                    end
                end

                %   check the cell is not too big, if requested
                if check_arg( maxSize )
                    if rp(i).Area > maxSize
                        valid = TOO_BIG;
                    end
                end

                %   check the cell is not too round, if requested
                if check_arg( aspectRatio )
                    if rp(i).MajorAxisLength / rp(i).MinorAxisLength <...
                            aspectRatio
                        valid = TOO_ROUND;
                    end
                end

                %   check the cell is solid enough, if requested
                if check_arg( solidity )
                    if rp(i).Solidity < solidity
                        valid = TOO_BENT;
                    end
                end

                %   check the cell is not too bright, if requested
                if check_arg( maxBright )
                    if rp(i).MeanIntensity > maxBright
                        valid = TOO_BRIGHT;
                    end
                end
                
                %   create a new cell object with all of the above
                self.cells{i} = ColiCell( rp(i).Image, rp(i).BoundingBox(1),...
                    rp(i).BoundingBox(2),...
                    rp(i).MeanIntensity - blackLevel, valid,...
                    self.imageSize, rpz(i).MeanIntensity, expansion );
               
            end
                        
        end
        
        function validCells = validCells( self, varargin )
            validCells = [];
            for i = 1:numel( self.cells )
                if self.cells{i}.validCell == 0
                    validCells = [ validCells, self.cells{i} ];
                end
            end
            if nargin == 2
                validCells = validCells(varargin{1});
            end
        end
        
        function getClusters( self, threshold, searchRadius, minSize,...
                maxSize )
            
            
            ballRadius = floor( 1000 / self.pixelSize );
            
            image = double( imread( self.clusterImage ) );
            blackLevel = getBlackLevel( self, 'cluster', 8 );
           
            cellBasalLevel = imopen( image, strel( 'ball', ballRadius,...
                ballRadius, 0 ) );
            
            statSet = statset( 'MaxIter', 25, 'TolFun', 1e-3,...
                'TolX', 1e-3 );
            
            minSize = minSize / self.pixelSize;
            maxSize = maxSize / self.pixelSize;
            searchRadius = floor( searchRadius / self.pixelSize );
            
            pixels = meshgrid( -searchRadius:searchRadius );
            
            xPix = reshape( pixels, 1, numel( pixels ) );
            yPix = reshape( pixels', 1, numel( pixels ) );
            
            validCells = self.validCells();
            
            for i = 1:numel( validCells )
                validCells(i).getClusters( image, cellBasalLevel,...
                    blackLevel, statSet, threshold, searchRadius,...
                    minSize, maxSize, xPix, yPix )
            end
            
        end
            
        function blackLevel = getBlackLevel( self, im, blackLevelCutoff )
            if strcmpi( im, 'cell' )
                image = imread( self.cellImage );
            elseif strcmpi( im, 'cluster' )
                image = imread( self.clusterImage );
            end
            image = image(:);
            blackLevel = mean( image(...
                image < ( 2^16 / blackLevelCutoff ) ) );
        end
        
        function displayOverlay( self )
            
            z = zeros( self.imageSize );
 
            
            figure(1);
            
            ceI = imread( self.cellImage );
            clI = imread( self.clusterImage );
            max(ceI(:))
            max(clI(:))
            im = cat( 3, ceI, clI, z ) ;
            max(im(:))
            size( im )
            
            imshow( cat( 3, ceI, clI, z ) );
            
        end
            
        
        function displayCells( self, im )
            
            load( 'errors.mat' );
            if strcmpi( im, 'cell' )
                image = imread( self.cellImage );
            elseif strcmpi( im, 'cluster' )
                image = imread( self.clusterImage );
            end
            
            cell_mask = zeros( self.imageSize );
            size_mask = zeros( self.imageSize );
            round_mask = zeros( self.imageSize );
            bright_mask = zeros( self.imageSize );
            bent_mask = zeros( self.imageSize );
            
            for i = 1:numel( self.cells )
                mask = self.cells{i}.expandedFullMask();
                if self.cells{i}.validCell == 0
                    cell_mask = cell_mask | mask;
                elseif self.cells{i}.validCell == TOO_BIG || self.cells{i}.validCell == TOO_SMALL
                    size_mask = size_mask | mask;
                elseif self.cells{i}.validCell == TOO_ROUND
                    round_mask = round_mask | mask;
                elseif self.cells{i}.validCell == TOO_BRIGHT
                    bright_mask = bright_mask | mask;
                elseif self.cells{i}.validCell == TOO_BENT
                    bent_mask = bent_mask | mask;
                end 
            end

            o = ones( self.imageSize );
            z = zeros( self.imageSize );

            red = cat( 3, o, z, z );
            green = cat( 3, z, o, z );
            blue = cat( 3, z, z, o );
            yellow = cat( 3, o, o, z );
            orange = cat( 3, o, 0.5 * o, z );
            

            figure(1);

            imshow(image, [ min( image(:) ), max( image(:) ) ] );
            hold on
            valid = imshow( green );
            set( valid, 'AlphaData', cell_mask * 192 );
            round = imshow( red );
            set( round, 'AlphaData', round_mask * 192 );
            bent = imshow( blue );
            set( bent, 'AlphaData', bent_mask * 192 );
            bright = imshow( yellow );
            set( bright, 'AlphaData', bright_mask * 128 );
            sizes = imshow( orange );
            set( sizes, 'AlphaData', size_mask * 128 );
            hold off;
            
        end
        
        function cc = clusterCount( self )
            validCells = self.validCells();
            cc = zeros( 1, numel( validCells ) );
            for i = 1:numel( validCells )
                cc( i ) = validCells(i).clusterCount();
            end
        end
        
        function varargout = clusterHistogram( self )
            h = histc( self.clusterCount(), 0:4 );
            varargout{1} = h;
            if nargout == 0
                bar( 0:4, h, 'histc' )
            end
        end
        
        function cf = clusterFraction( self )
            validCells = self.validCells();
            cf = [];
            for i = 1:numel( validCells )
                cf = [ cf, validCells(i).clusterFraction ];
            end
        end
        
            
    end
end

function ok = check_arg( arg )

ok = isa( arg, 'double' ) && numel( arg ) == 1 && arg > 0;

end
