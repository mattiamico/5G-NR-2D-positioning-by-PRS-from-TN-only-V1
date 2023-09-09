function gNBPos = getgNB2DPositions(numgNBs, r1, r2, central_position, angularlyEquispaced)
%   GNBPOS = getgNBPositions(NUMGNBS) returns a cell array of random
%   positions for gNBs, given the number of gNBs NUMGNBS.     
    arguments
        numgNBs (1,1) {mustBeInteger} %  = 4
        r1 (1,1) double % 3000 [m]
        r2 (1,1) double % 4000 [m]
        central_position (1,2) double = [0 0]
        angularlyEquispaced (1,1) logical = 0 
    end

    gNBPos = cell(1,numgNBs);
    for gNBIDx = 1:numgNBs
        
        if angularlyEquispaced
            % position gNBs at equispaced radial positions if chosen 
            if mod(numgNBs,2)==0 
                phis = pi/(numgNBs) + (0:numgNBs)*2*pi/numgNBs;
            else 
                phis = pi/(2*numgNBs) + (0:numgNBs)*2*pi/numgNBs;
            end
            phi = phis(gNBIDx);
        else
            % Position gNB randomly within gNBIDx*2*pi/numgNBs radian sector/slice      
            phi = gNBIDx*2*pi/numgNBs + rand(1,1)*2*pi/(2*numgNBs) - 2*pi/(2*numgNBs);
        end

        % Position gNB randomly between r1 and r2 from UE aka r falls between [r1,r2]
        r = randi([0,r2-r1],1,1) + r1;
        s = rng;
    
        % Convert polar coordinates to Cartesian coordinates
        [x,y] = pol2cart(phi,r);
        x_shifted = x+central_position(1);
        y_shifted = y+central_position(2);

        %gNBPos{gNBIDx} = [x,y];
        gNBPos{gNBIDx} = [x_shifted,y_shifted]; 

    end

end









