function plotWaveforms(waveforms,carrier,title_prefix)

    numgNBs = numel(waveforms);
    colors = getColors(numgNBs);
    figure;
    hold on;
    
    % Plot correlation for each gNodeB
    t = 0:length(waveforms{1})-1;
    t = t';
    legendstr = cell(1,2*numgNBs);

    for gNBIdx = 1:numgNBs
        plot(t,abs(waveforms{gNBIdx}), ...
            'Color',colors(gNBIdx,:),'LineWidth',2);

        % legendstr{gNBIdx} = sprintf('gNB%d (NCellID = %d)', ...
        %     gNBIdx,carrier(gNBIdx).NCellID);
    end
    
    %legend(legendstr);
    xlabel('Time [s]')
    ylabel('Signal abs value');
    title(title_prefix + ' waveforms');
end








