


function plotPositionsAndHyperbolaCurves(gNBPos,UEPos,detgNBNums,curveX,curveY,gNBNums,estPos)
%   plotPositionsAndHyperbolaCurves(GNBPOS,UEPOS,DETGNBNUMS,CURVEX,CURVEY,GNBNUMS,ESTPOS)
%   plots gNB, UE positions, and hyperbola curves by considering these
%   inputs:
%   GNBPOS     - Positions of gNBs
%   UEPOS      - Position of UE
%   DETGNBNUMS - Detected gNB numbers
%   CURVEX     - Cell array having the x-coordinates of hyperbola curves
%   CURVEY     - Cell array having the y-coordinates of hyperbola curves
%   GNBNUMS    - Cell array of gNB numbers corresponding to all hyperbolas
%   ESTPOS     - Estimated UE position


    plotgNBAndUE_2DPositions(gNBPos,UEPos,detgNBNums, [], [], [], "", false, []); % edited wrt original

    for curveIdx = 1:numel(curveX)
        curves(curveIdx) = plot(curveX{1,curveIdx},curveY{1,curveIdx}, ...
               '--','LineWidth',1,'Color','k');
    end
    % Add estimated UE position to figure
    plot(estPos(1),estPos(2),'+','MarkerSize',10, ...
        'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','LineWidth',2);

    % Add legend
    gNBLegends = cellstr(repmat("gNB",1,numel(detgNBNums)) + ...
        detgNBNums + [" (reference gNB)" repmat("",1,numel(detgNBNums)-1)]);
    legend([gNBLegends {'Actual UE Position'} repmat({''},1,numel(curveX)) {'Estimated UE Position'}]);

    for curveIdx = 1:numel(curves)
        % Create a data tip and add a row to the data tip to display
        % gNBs information for each hyperbola curve
        dt = datatip(curves(curveIdx));
        row = dataTipTextRow("Hyperbola of gNB" + gNBNums{curveIdx}(1) + ...
            " and gNB" + gNBNums{curveIdx}(2),'');
        curves(curveIdx).DataTipTemplate.DataTipRows = [row; curves(curveIdx).DataTipTemplate.DataTipRows];
        % Delete the data tip that is added above
        delete(dt);
    end

end