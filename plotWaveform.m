function plotWaveform(waveform, sr, title_prefix, carrier, gNBIdx, numgNBs)

    colors = getColors(numgNBs);
    figure;
    hold on;
    
    t = (0:length(waveform) -1)/sr;    
    plot(t,abs(waveform), ...
        'Color',colors(gNBIdx,:),'LineWidth',2);

    legendstr{gNBIdx} = sprintf('gNB%d (NCellID = %d)', ...
         gNBIdx,carrier(gNBIdx).NCellID);
    
    legend(legendstr);
    xlabel('Time [s]')
    ylabel('Signal Absolute value');
    title(title_prefix + 'gNB' + gNBIdx);

    
    % TODO: if the signal is complex: plot also the real and imag parts as
    % 2 other subplots ( aka 2 subplots rows )

end

































