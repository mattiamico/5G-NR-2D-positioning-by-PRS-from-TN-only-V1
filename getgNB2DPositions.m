function gNBPos = getgNB2DPositions(numgNBs, r1, r2)

%   GNBPOS = getgNBPositions(NUMGNBS) returns a cell array of random
%   positions for gNBs, given the number of gNBs NUMGNBS.
    
    gNBPos = cell(1,numgNBs);
    for gNBIDx = 1:numgNBs
        % Position gNB randomly within gNBIDx*2*pi/numgNBs radian sector/slice
        phi = gNBIDx*2*pi/numgNBs + rand(1,1)*2*pi/(2*numgNBs) - 2*pi/(2*numgNBs);

        % Position gNB randomly between r1+(gNBIDx*r2/numgNBs) and r2+(gNBIDx*r2/numgNBs) from UE
        r = randi([0,r2-r1],1,1) + r1 + (gNBIDx*r2/numgNBs);
        s = rng;
    
        % Convert polar coordinates to Cartesian coordinates
        [x,y] = pol2cart(phi,r);
        gNBPos{gNBIDx} = [x,y];
    end

end










