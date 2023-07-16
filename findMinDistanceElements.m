function [firstCurvePoints,secondCurvePoints] = findMinDistanceElements(xA,yA,xB,yB)
%   [FIRSTCURVEPOINTS,SECONDCURVEPOINTS] = findMinDistanceElements(XA,YA,XB,YB)
%   returns the closest points between the given hyperbolic surfaces.

    distAB1 = zeros(numel(xA),numel(xB));
    for idx1 = 1:numel(xA)
        distAB1(idx1,:) = sqrt((xB-xA(idx1)).^2 + (yB-yA(idx1)).^2);
    end
    [~,rows] = min(distAB1,[],'omitnan');
    [~,col] = min(min(distAB1,[],'omitnan'));

    % Extract the indices of the points on two curves which are up to 5
    % meters apart. This is to identify the multiple intersections between
    % two hyperbola curves.
    [allRows,allCols] = find(distAB1 <= distAB1(rows(col),col)+5);
    diffVec = abs(allRows - allRows(1));
    index = find(diffVec > (min(diffVec) + max(diffVec))/2,1);
    allRows = [allRows(1);allRows(index)];
    allCols = [allCols(1);allCols(index)];

    for idx = 1:numel(allRows)
        firstCurveIndices = allRows(idx);
        secondCurveIndices = allCols(idx);
        x1a = xA(firstCurveIndices);
        y1a = yA(firstCurveIndices);
        x2a = xB(secondCurveIndices);
        y2a = yB(secondCurveIndices);

        % Use subsequent points to create line for linearization
        if firstCurveIndices == numel(rows)
            x1b = xA(firstCurveIndices-1);
            y1b = yA(firstCurveIndices-1);
        else
            x1b = xA(firstCurveIndices+1);
            y1b = yA(firstCurveIndices+1);
        end
        if secondCurveIndices == numel(rows)
            x2b = xB(secondCurveIndices-1);
            y2b = yB(secondCurveIndices-1);
        else
            x2b = xB(secondCurveIndices+1);
            y2b = yB(secondCurveIndices+1);
        end
        firstCurvePoints{idx} = [x1a y1a;x1b y1b]; %#ok<*AGROW> 
        secondCurvePoints{idx} = [x2a y2a;x2b y2b];
    end
end