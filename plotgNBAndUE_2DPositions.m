function plotgNBAndUE_2DPositions(gNBPos,UEPos,gNBNums,r1,r2,central_position, title_string, savePlot,outputDirectory)
%   plotgNBAndUE_2DPositions(GNBPOS,UEPOS,GNBNUMS) plots gNB and UE 2D positions.
    
    numgNBs = numel(gNBNums);
    colors = getColors(numel(gNBPos));
    
    f = figure;
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
        % Plot r1 circumference and r2 circumference
        % center = [ 0 0; 0 0];
        center = [central_position(1:2) ; central_position(1:2)];
        viscircles(center, [r1 r2],'Color','r','LineStyle','--');
    end
   
    axis equal;
    legend([lengendstr 'UE'])
    title(title_string,'FontSize',11,'Interpreter','latex')
    xlabel('X Position [m]')
    ylabel('Y Position [m]')

    grid on;
    
    if savePlot
        %set(f, 'PaperPositionMode', 'auto')
        set(gca,'YMinorTick','on','YMinorGrid','on');
        set(gca,'XMinorTick','on','XMinorGrid','on','XTickMode','manual');
        saveas(f, fullfile(outputDirectory,'gNBAndUE_2DPositions'),'png')
    end 
end


