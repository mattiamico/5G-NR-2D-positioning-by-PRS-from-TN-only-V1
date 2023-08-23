function gNBPos = getgNBPositions(numgNBs)
%   GNBPOS = getgNBPositions(NUMGNBS) returns a cell array of random
%   positions for gNBs, given the number of gNBs NUMGNBS (AS IN THE
%   ORIGINAL MATLAB TUTORIAL)

%   NB: gNBs are displaced circularly around the cartesian plane center and
%       not around the UE position ==> ISSUE??

    gNBPos = cell(1,numgNBs);
    for gNBIdx = 1:numgNBs
        % Position gNB randomly within gNBNum*2*pi/numgNBs radian sector
        phi = gNBIdx*2*pi/numgNBs + rand(1,1)*2*pi/(2*numgNBs) - 2*pi/(2*numgNBs);

        % Position gNB randomly between 4000 + (gNBNum*5000/numgNBs) and 5000 +
        % (gNBNum*5000/numgNBs) from UE
        r = randi([0,1000],1,1) + 4000 + (gNBIdx*5000/numgNBs);

        % Convert polar coordinates to Cartesian coordinates
        [x,y] = pol2cart(phi,r);
        gNBPos{gNBIdx} = [x,y];
    end
end