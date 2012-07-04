
classdef ColiCell < handle

    properties
       
        validCell
        cellMask
        topLeftX
        topLeftY
        parentSize
        area
        length
        tilt
        cheYMeanIntensity
        cheZMeanIntensity
        clusters
        expansion
        
    end
    
    methods
        
        function self = ColiCell( cellMask, tlx, tly, mIntensity,...
                valid, sz, mIntensityZ, expansion )
            self.cellMask = zeros( size( cellMask ) +...
                ( 2 * [ expansion expansion ] ) );
            self.cellMask( expansion + 1:end - expansion,...
                expansion + 1:end-expansion ) = cellMask;
            
            self.validCell = valid;
            self.topLeftX = ceil( tlx ) - expansion;
            self.topLeftY = ceil( tly ) - expansion;
            self.cheYMeanIntensity = mIntensity;
            self.cheZMeanIntensity = mIntensityZ;
            self.parentSize = sz;
            self.expansion = expansion;

            
            if self.validCell == 0
                rp = regionprops( cellMask, 'Area',...
                    'MajorAxisLength', 'Orientation' );
                self.area = rp.Area;
                self.length = rp.MajorAxisLength;
                self.tilt = rp.Orientation;
            end
        end
        
        function mask = expandedMask( self )
            mask = bwmorph( self.cellMask, 'dilate',...
                floor( self.expansion / 1.5 ) );
        end
        
        function mask = fullMask( self )
            mask = zeros( self.parentSize );
            mask( self.topLeftY:self.topLeftY + size( self.cellMask, 1 ) - 1,...
                self.topLeftX:self.topLeftX + size( self.cellMask, 2 ) - 1 )...
                = self.cellMask;
        end
        
        function mask = expandedFullMask( self )
            mask = zeros( self.parentSize );
            mask( self.topLeftY:self.topLeftY + size( self.cellMask, 1 ) - 1,...
                self.topLeftX:self.topLeftX + size( self.cellMask, 2 ) - 1 )...
                = self.expandedMask();
        end
        
        function image = subImage( self, image )
            image = image( self.topLeftY:self.topLeftY +...
                size( self.cellMask, 1 ) - 1,...
                self.topLeftX:self.topLeftX +...
                size( self.cellMask, 2 ) - 1 );
        end
        
        
        function getClusters( self, image, cellBasalLevel, blackLevel,...
                statSet, threshold, searchRadius, minSize, maxSize,...
                xPix, yPix )
                        
            maskedImage = image .* self.expandedFullMask();
            figure(1);
            imshow( maskedImage, [min(maskedImage(:)), max(maskedImage(:))]);
            self.clusters = {};
            
            count = 0;
            
            while(1)
                ok = true;
                count = count + 1;
                                
                [ maxPixelValue, maxPixelIndex ] = max( maskedImage(:) );
                [ maxPixelY, maxPixelX ] = ind2sub( size(...
                    maskedImage ), maxPixelIndex );
                                
                if ( maxPixelValue < threshold * self.cheZMeanIntensity )
                    break
                end
                               
                
                clusterBoxXMin = maxPixelX - searchRadius;
                clusterBoxXMax = maxPixelX + searchRadius;
                
                clusterBoxYMin = maxPixelY - searchRadius;
                clusterBoxYMax = maxPixelY + searchRadius;
                
                [ clusterBoxY, clusterBoxX ] = meshgrid(...
                    clusterBoxYMin:clusterBoxYMax,...
                    clusterBoxXMin:clusterBoxXMax );
                
                clusterImage = image(...
                    clusterBoxYMin:clusterBoxYMax,...
                    clusterBoxXMin:clusterBoxXMax );
                
                clusterBasalLevel = cellBasalLevel(...
                    clusterBoxYMin:clusterBoxYMax,...
                    clusterBoxXMin:clusterBoxXMax );
                
                sumMass = sum( clusterImage(:) );
                
                sumXMass = sum( sum( clusterBoxX .* clusterImage, 1 ), 2 );
                sumYMass = sum( sum( clusterBoxY .* clusterImage, 1 ), 2 );
                
                clusterCoMX = sumXMass / sumMass;
                clusterCoMY = sumYMass / sumMass;
                
                % CENTROID TO BE CENTERED ON CLUSTER
                % centroid is given relative to maxPixelIndex always
                
                centroidX = clusterCoMX - maxPixelX;
                centroidY = clusterCoMY - maxPixelY;
                
                beta = [ centroidX, centroidY, maxPixelValue * 1.5, 1 ];
                
                xyInput = zeros( 1, numel( clusterImage ) );
                zOutput = double( reshape( clusterImage -...
                    clusterBasalLevel - blackLevel,...
                    1, numel( clusterImage ) ) );
                
                beta = nlinfit( xyInput, zOutput, @( beta, ~ )...
                    beta(3) *  exp( -2 * ( ( xPix - beta(1) ) .^ 2 +...
                    ( yPix - beta(2) ) .^ 2 ) /...
                    ( beta(4) ) ^ 2 ), beta, statSet );
                
                self.clusters{ count } = ColiCluster(...
                    beta(1) + maxPixelX, beta(2) + maxPixelY,...
                    beta(3), beta(4) );
                
                
                maskedImage(...
                    clusterBoxYMin:clusterBoxYMax,...
                    clusterBoxXMin:clusterBoxXMax ) = 0;
                
     
                
                
            end
            
        end
        
        function cc = clusterCount( self )
            cc = numel( self.clusters );
        end
        
    end
    
end
