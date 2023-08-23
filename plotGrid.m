function plotGrid(prsGrid,dataGrid)
%   plotGrid(PRSGRID,DATAGRID) plots the carrier grid from all gNBs.

    numgNBs = numel(prsGrid);
    figure()
    mymap = [1 1 1; ...       % White color for background
             0.8 0.8 0.8; ... % Gray color for PDSCH data from all gNBs
             getColors(numgNBs)];
    chpval = 3:numgNBs+2;

    gridWithPRS = zeros(size(prsGrid{1}));
    gridWithData = zeros(size(dataGrid{1}));
    names = cell(1,numgNBs);
    for gNBIdx = 1:numgNBs
        gridWithPRS = gridWithPRS + abs(prsGrid{gNBIdx})*chpval(gNBIdx);
        gridWithData = gridWithData + abs(dataGrid{gNBIdx});
        names{gNBIdx} = ['PRS from gNB' num2str(gNBIdx)];
    end
    names = [{'Data from all gNBs'} names];
    % Replace all zero values of gridWithPRS with proper scaling (value of 1)
    % for white background
    gridWithPRS(gridWithPRS == 0) = 1;
    gridWithData(gridWithData ~=0) = 1;

    % Plot grid
    image(gridWithPRS+gridWithData); % In the resultant grid, 1 represents
                                     % the white background, 2 (constitute
                                     % of 1s from gridWithPRS and
                                     % gridWithData) represents PDSCH, and
                                     % so on
    % Apply colormap
    colormap(mymap);
    axis xy;

    % Generate lines
    L = line(ones(numgNBs+1),ones(numgNBs+1),'LineWidth',8);
    % Index colormap and associate selected colors with lines
    set(L,{'color'},mat2cell(mymap(2:numgNBs+2,:),ones(1,numgNBs+1),3)); % Set the colors

    % Create legend
    legend(names{:});

    % Add title
    title('Carrier Grid Containing PRS and PDSCH from Multiple gNBs');

    % Add labels to x-axis and y-axis
    xlabel('OFDM Symbols');
    ylabel('Subcarriers');
end
