function plotgNBAndUE_2DPositions(gNBPos,UEPos,gNBNums,r1,r2)
%   plotgNBAndUE_2DPositions(GNBPOS,UEPOS,GNBNUMS) plots gNB and UE 2D positions.
    
    numgNBs = numel(gNBNums);
    colors = getColors(numel(gNBPos));
    
    figure;
    hold on;
    lengendstr = cell(1,numgNBs);
    
    % Plot position of each gNB
    for gNBIdx = 1:numgNBs
        plot(gNBPos{gNBNums(gNBIdx)}(1), gNBPos{gNBNums(gNBIdx)}(2), ...
             'Color',colors(gNBNums(gNBIdx),:),'Marker','^', ...
             'MarkerSize',11,'LineWidth',2,'MarkerFaceColor',colors(gNBNums(gNBIdx),:));
        lengendstr{gNBIdx} = sprintf('gNB%d',gNBIdx);
    end 

    % Plot UE position
    plot(UEPos(1),UEPos(2), 'ko', 'MarkerSize', 10, 'LineWidth',2, 'MarkerFaceColor','k');

    if ~(isempty(r1) && isempty(r2))
        % TODO: Plot r1 circumference and r2 circumference
        center = [ 0 0; 0 0];
        viscircles(center, [r1 r2],'Color','r','LineStyle','--');
    end
   
    axis equal;
    legend([lengendstr 'UE'])
    xlabel('X Position [m]')
    ylabel('Y Position [m]')

end


