function gNBPos = getgNBPositions(numgNBs)

%   GNBPOS = getgNBPositions(NUMGNBS) returns a cell array of random
%   positions for gNBs, given the number of gNBs NUMGNBS.
    
    gNBPos = cell(1,numgNBs);
    for gNBIDx = 1:numgNBs
        % Position gNB randomly within gNBIDx*2*pi/numgNBs radian sector/slice
        phi = gNBIDx*2*pi/numgNBs + rand(1,1)*2*pi/(2*numgNBs) - 2*pi/(2*numgNBs);

        % Position gNB randomly between 4000 + (gNBIDx*5000/numgNBs) and 5000 +
        % (gNBIDx*5000/numgNBs) from UE
        r = randi([0,1000],1,1) + 4000 + (gNBIDx*5000/numgNBs);
    
        % Convert polar coordinates to Cartesian coordinates
        [x,y] = pol2cart(phi,r);
        gNBPos{gNBIDx} = [x,y];
    end





end
