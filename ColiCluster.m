
classdef ColiCluster < handle
    
    properties
        
        centroidX
        centroidY
        radius
        height
        
    end
    
    methods
        
        function self = ColiCluster( cX, cY, h, r )
            
            self.centroidX = cX;
            self.centroidY = cY;
            self.radius = r;
            self.height = h;
            
        end
        
        function volume = volume( self )
            
            volume = self.radius ^ 2 * self.height;
            
        end
        
        function fraction = fraction( self, meanIntensity, area )
            fraction = self.volume() / ( meanIntensity * area );
        end
        
    end
    
end