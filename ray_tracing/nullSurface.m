classdef nullSurface < opticalElement_Surface
% Null surface used for initialization purposes in class
% opticalElement_Surface

    methods
        function goToSurf(~,~)
            error('opticalElement_Surface:nullSurface','A nullSurface may not be used in calculations. Look for a problem with initialization.')
        end
        function createSurfPatch(~,~,~)
            error('opticalElement_Surface:nullSurface','A nullSurface may not be used in calculations. Look for a problem with initialization.')
        end
        function flipSurf(~,~,~)
            error('opticalElement_Surface:nullSurface','A nullSurface may not be used in calculations. Look for a problem with initialization.')
        end
    end
    
end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
