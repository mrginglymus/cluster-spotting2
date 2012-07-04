
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
        
    end
    
    methods
        
        function self = ColiCell( cellMask, tlx, tly, mIntensity,...
                valid, sz, mIntensityZ )
            self.cellMask = cellMask;
            self.validCell = valid;
            self.cellMask = cellMask;
            self.topLeftX = ceil( tlx );
            self.topLeftY = ceil( tly );
            self.cheYMeanIntensity = mIntensity;
            self.cheZMeanIntensity = mIntensityZ;
            self.parentSize = sz;

            
            if self.validCell == 0
                rp = regionprops( cellMask, 'Area',...
                    'MajorAxisLength', 'Orientation' );
                self.area = rp.Area;
                self.length = rp.MajorAxisLength;
                self.tilt = rp.Orientation;
            end
        end
        
        function mask = fullMask( self )
            mask = zeros( self.parentSize );
            mask( self.topLeftY:self.topLeftY + size( self.cellMask, 1 ) - 1,...
                self.topLeftX:self.topLeftX + size( self.cellMask, 2 ) - 1 )...
                = self.cellMask;
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
            
            maskedImage = image .* self.fullMask();
            figure(1);
            imshow( maskedImage, [min(maskedImage(:)), max(maskedImage(:))]);
            self.clusters = [];
            
            
            while(1)
                ok = true;
                                
                [ maxPixelValue, maxPixelIndex ] = max( maskedImage(:) );
                [ maxPixelY, maxPixelX ] = ind2sub( size(...
                    maskedImage ), maxPixelIndex );
                
                disp( maxPixelValue )
                disp( self.cheZMeanIntensity )
                
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
                
                % CENTROID TO BE CENTERED ON CLUSTER (still in cell space)
                % centroid is given relative to maxPixelIndex always
                
                centroidX = clusterCoMX - maxPixelX;
                centroidY = clusterCoMY - maxPixelY;
                
                beta = [ centroidX, centroidY, maxPixelValue * 1.5, 1 ];
                
                xyInput = zeros( 1, numel( clusterImage ) );
                zOutput = double( reshape( clusterImage -...
                    clusterBasalLevel, 1, numel( clusterImage ) ) );
                figure(5000)
                mesh(clusterImage)
                %{
                beta = nlinfit( xyInput, zOutput, @( beta, ~ )...
                    beta(3) *  exp( -2 * ( ( xPix - beta(1) ) .^ 2 +...
                    ( yPix - beta(2) ) .^ 2 ) /...
                    ( beta(4) ) ^ 2 ), beta, statSet );
                %}
                
                maskedImage(...
                    clusterBoxYMin:clusterBoxYMax,...
                    clusterBoxXMin:clusterBoxXMax ) = 0;
                
                figure(2);
                imshow(maskedImage);
                
                pause
                
            end
        end
        
    end
    
end


function z = gaussianMesh( beta, ~ )
global xPix yPix;

z = beta(3)*exp(-2*((xPix-beta(1)).^2+(yPix-beta(2)).^2)/(beta(4))^2);
z = reshape( z, 1, numel(z) );
end