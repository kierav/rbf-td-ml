function xt = create_nodes_01(rFactor,gridtype,Tf)
% create_nodes_01 generates nodes in the unit square of the x,t-plane
% --- Input parameters ---
%   rFactor     Desired distance between nodes
%   gridtype    String indicating type of nodeset (quasiunif, hex, cart,
%                   cart2, vardensity)
%   Tf          Final time
% --- Output parameter ---
%   xt          All nodes in the x,t-plane, sorted for increasing t-values

if nargin < 2
    gridtype = 'quasiunif';
    Tf = 1;
end

% Create the scattered node set in the x,t-plane
switch gridtype
    case 'quasiunif'
        dx = radius_2d([0,0],rFactor);
        ninit = round(1./dx*10);
        dotmax = round((1./dx)^2*2);
        xt = node_drop_2d([0 1 0 Tf],ninit,dotmax,@radius_2d,rFactor);
        
        % Create boundary nodes
        corners = [0 0;1 0;1 Tf;0 Tf];        % The four corners
        bdyr = discretize_bdy (corners,@radius_2d,rFactor);
        bdy = bdyr(:,1:2);
        xt2  = repel(xt,bdy,corners,@radius_2d,rFactor); % Repel by all boundaries
    case 'hex'
        % lay down hex grid
        dx = rFactor;
        x1 = 0:dx:1;
        x2 = (0+dx/2):dx:1;
        y = 0:dx*sqrt(3)/2:Tf;

        x = [x1';x2'];
        X = repmat(x,floor(length(y)/2),1);
        if mod(length(y),2) ~= 0
            X = [X;x1'];
        end
        Y = [];
        for i = 1:length(y)
            if mod(i,2) == 0
                Y = [Y;y(i)*ones(length(x1),1)];
            else
                Y = [Y;y(i)*ones(length(x2),1)];
            end
        end  
        Ibdy = X == 0 | X == 1 | Y == 0 | Y == Tf | Y == max(y);
        xt2 = [X(~Ibdy) Y(~Ibdy)];
        bdy = [X(Ibdy) Y(Ibdy)];
    case 'cart'
        % lay down Cartesian grid
        dx = rFactor;
        x = 0:dx:1;
        y = 0:dx:Tf;
        [X,Y] = meshgrid(x,y);
        X = X(:); Y = Y(:); 
        Ibdy = X == 0 | X == 1 | Y == 0 | Y == Tf;
        xt2 = [X(~Ibdy) Y(~Ibdy)];
        bdy = [X(Ibdy) Y(Ibdy)];
    case 'cart2'
        % lay down Cartesian grid with stretched t coordinate
        dx = rFactor;
        dt = rFactor^2*15;
        x = 0:dx:1;
        y = 0:dt:1;
        [X,Y] = meshgrid(x,y);
        X = X(:); Y = Y(:); 
        Ibdy = X == 0 | X == 1 | Y == 0 | Y == 1;
        xt2 = [X(~Ibdy) Y(~Ibdy)];
        bdy = [X(Ibdy) Y(Ibdy)];
    case 'vardensity'
        ninit = 1e6;
        dotmax = 1e7;
        xt = node_drop_2d([0 1 0 1],ninit,dotmax,@radius_1dTD,rFactor);
        
        % Create boundary nodes
        corners = [0 0;1 0;1 1;0 1];        % The four corners
        bdyr = discretize_bdy(corners,@radius_1dTD,rFactor);
        bdy = bdyr(:,1:2);
        xt2  = repel(xt,bdy,corners,@radius_1dTD,rFactor); % Repel by all boundaries
end

xt = xt2;   % Remove last column with radius information

[ni,~] = size(xt);      % Get total number of interior nodes
[nb,~] = size(bdy);     % Get total number of boundary nodes

xt = [bdy;xt];          % Make xy contain both boundary and interior nodes
xt = sortrows(xt,2);    % Sort the rows from lowest y (i.e. t) to highest