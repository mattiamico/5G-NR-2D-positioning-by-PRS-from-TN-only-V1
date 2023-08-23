function validateNumLayers(pdsch)
%   validateNumLayers(PDSCH) validates the number of transmission layers,
%   given the array of physical downlink shared channel configuration
%   objects PDSCH.
    
    numLayers = [pdsch(:).NumLayers];
    if ~all(numLayers == 1)
        error('nr5g:invalidNLayers',['The number of transmission layers ' ...
            'configured for the data transmission must be 1.']);
    end
end