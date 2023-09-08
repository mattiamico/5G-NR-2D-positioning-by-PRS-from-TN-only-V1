function gNBPos = getgNB2DPositions(numgNBs, r1, r2)

%   GNBPOS = getgNBPositions(NUMGNBS) returns a cell array of random
%   positions for gNBs, given the number of gNBs NUMGNBS.
    
%   NB: gNBs are displaced circularly around the cartesian plane center and
%       not around the UE position ==> ISSUE??

    gNBPos = cell(1,numgNBs);
    for gNBIDx = 1:numgNBs
        
        % Position gNB randomly within gNBIDx*2*pi/numgNBs radian sector/slice % FIXME: position at equispaced positions instead ?? 
        phi = gNBIDx*2*pi/numgNBs + rand(1,1)*2*pi/(2*numgNBs) - 2*pi/(2*numgNBs);

        % Position gNB randomly between r1 and r2 from UE aka r falls between [r1,r2]
        r = randi([0,r2-r1],1,1) + r1;
        s = rng;
    
        % Convert polar coordinates to Cartesian coordinates
        [x,y] = pol2cart(phi,r);
        gNBPos{gNBIDx} = [x,y];
    end

end










