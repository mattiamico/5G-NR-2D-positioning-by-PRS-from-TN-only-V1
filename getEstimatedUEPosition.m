function estPos = getEstimatedUEPosition(xCell,yCell)
%   ESTPOS = getEstimatedUEPosition(XCELL,YCELL) returns the x- and y-
%   coordinates of the UE position, given the hyperbolas XCELL and YCELL.

    % Steps involved in this computation:
    % 1. Find closest points between hyperbolic surfaces
    % 2. Linearize surfaces around the points closest to intersection to
    % interpolate actual intersection location

    numCurves = numel(xCell);
    % Make all vectors have equal length
    maxLen = max(cellfun(@length,xCell));
    for curveIdx = 1:numCurves
        xCell{curveIdx} = [xCell{curveIdx} inf(1,maxLen-length(xCell{curveIdx}))];
        yCell{curveIdx} = [yCell{curveIdx} inf(1,maxLen-length(yCell{curveIdx}))];
    end
    tempIdx = 1;
    for idx1 = 1:numCurves-1
        for idx2 = (idx1+1):numCurves
            [firstCurve,secondCurve] = findMinDistanceElements(xCell{idx1},yCell{idx1},xCell{idx2},yCell{idx2});
            for idx3 = 1:numel(firstCurve)
                [x1a,y1a,x1b,y1b] = deal(firstCurve{idx3}(1,1),firstCurve{idx3}(1,2), ...
                                         firstCurve{idx3}(2,1),firstCurve{idx3}(2,2));
                [x2a,y2a,x2b,y2b] = deal(secondCurve{idx3}(1,1),secondCurve{idx3}(1,2), ...
                                         secondCurve{idx3}(2,1),secondCurve{idx3}(2,2));
                a1 = (y1b-y1a)/(x1b-x1a);
                b1 = y1a - a1*x1a;

                a2 = (y2b-y2a)/(x2b-x2a);
                b2 = y2a - a2*x2a;

                xC(tempIdx) = (b2-b1)/(a1-a2);
                yC(tempIdx) = a1*xC(tempIdx) + b1;
                tempIdx = tempIdx+1;
            end
        end
    end
    estPos = [mean(xC) mean(yC)];
end
