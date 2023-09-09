


function [x,y,delta,phi,r, hk] = getRSTDCurve_improved(gNB1,gNB2,rstd)

%   [X,Y] = getRSTDCurve(GNB1,GNB2,RSTD) returns the x- and y-coordinates
%   of a hyperbola equation corresponding to an RSTD value [m], given the pair
%   of gNB positions gNB1 and gNB2.

    %% TD PLOT FOR DEBUG & DEV
    % Calculate vector between two gNBs
    delta = gNB1 - gNB2;

    
    % Express distance vector in polar form
    [phi,r] = cart2pol(delta(1),delta(2));  % r = 2c

    r_ref = sqrt(delta(1)^2 + delta(2)^2);


    rd = (r+rstd)/2;        % distance of segment between Vertex and Focus 

    % Compute the hyperbola parameters

    a = (r/2)-rd;           % Vertex      a = c - rd 

    c = r/2;                % Focus       c = a + rd

    b = sqrt(c^2-a^2);      % Co-vertex


    %% TD PLOT FOR DEBUG & DEV
    hk = (gNB1 + gNB2)/2;   % Midpoint vector between the two gNBs ( ie center of hyperbola aka vettore di shift f0 di wikipedia parametric def  )      
    
    mu = -2:1e-3:2; % parameter defining hyperbola extension/truncation 
    %mu = -3:1e-3:3; % parameter defining hyperbola extension/truncation 


    % Get x- and y- coordinates of (parametric) hyperbola equation

    x = (a*cosh(mu)*cos(phi)-b*sinh(mu)*sin(phi)) + hk(1);
    y = (a*cosh(mu)*sin(phi)+b*sinh(mu)*cos(phi)) + hk(2);

    % NB: right branch of the unit hyperbola (+- cosh(mu),sinh(mu))
    %x = cosh(mu);
    %y = sinh(mu);
    
    % NB: left branch of the non-unit hyperbola (+- a*cosh(mu),b*sinh(mu))
    %x = a*cosh(mu) + hk(1);
    %y = b*sinh(mu) + hk(2);
    


    %% DONE: DERIVARE VETTORI colonna f1 e f2 ( aka matrice di rotazione A di wikipedia parametric def )
    %% DONE: COMPUTARE x e y da calcoli matriciali
    
    non_unit_hyperbola = [a*cosh(mu) ; b*sinh(mu)];
    A = [cos(phi) -sin(phi) ; sin(phi) cos(phi)];

    hyp_ref = A*non_unit_hyperbola + hk';


    %x_ref = [ cos(pi) -sin(phi) ]*non_unit_hyp +  hk(1);
    eps_x = max(abs(x - hyp_ref(1,:))) 
    eps_y = max(abs(y - hyp_ref(2,:))) 

    disp("val ok?")

end
 






















