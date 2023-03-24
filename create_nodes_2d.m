function xt = create_nodes_2d(rFactor,gridtype,Tf)
% create_nodes_2d generates nodes in x,y,t space
% --- Input parameters ---
%   rFactor     Desired distance between nodes
%   gridtype    String indicating type of nodeset (quasiunif, cart, lake, vardensity)
%   Tf          Final time
% --- Output parameter ---
%   xt          All nodes in x,y,t-space, sorted for increasing t-values


if nargin < 3
    Tf = 1;
end
if nargin < 2
    gridtype = 'quasiunif';
end

% Create the scattered node set in x,y,t-space

switch gridtype
    case 'quasiunif'
        dx = radius_3d_unif([0,0,0],rFactor);
        ninit = round(1./dx*10);
        dotmax = round((1./dx)^2*(Tf./dx)*2);
        xt = node_drop_3d([0 1 0 1 0 Tf],[ninit,ninit],dotmax,@radius_3d_unif,rFactor);
        
        % Create boundary nodes
        nbdy = floor(1/dx);
        if 1-nbdy*dx > 0.2*dx
            nbdy = nbdy+1;
        end
        dxbdy = 1/nbdy;
        xbdy = (0:dxbdy:1)';
        ntbdy = round(Tf/dx);
        tbdy = linspace(0,Tf,ntbdy+1)';
        [XX,YY,TT] = meshgrid(xbdy,xbdy,tbdy);
        Ibdy = XX==0 | XX==1 | YY==0 | YY==1 | TT==Tf;
        XX = XX(Ibdy);
        YY = YY(Ibdy);
        TT = TT(Ibdy);
        bdy = [XX(:),YY(:),TT(:)];
        xt2  = repel_3d(xt,bdy,@radius_3d_unif,rFactor); % Repel by all boundaries
        xt2 = xt2(xt2(:,1)<1 & xt2(:,2)<1 & xt2(:,3)<Tf & xt2(:,1)>0 & xt2(:,2)>0 & xt2(:,3)>0,:);
    case 'vardensity'
        dx = rFactor;
        ninit = round(1./dx*10);
        dotmax = round((1./dx)^2*(Tf./dx)*3);
        xt = node_drop_3d([0 1 0 1 0 Tf],[ninit,ninit],dotmax,@radius_2dTD,rFactor);
        
        % Create boundary nodes
        nbdy = floor(1/dx);
        if 1-nbdy*dx > 0.2*dx
            nbdy = nbdy+1;
        end
        dxbdy = 1/nbdy;
        xbdy = (0:dxbdy:1)';
        ntbdy = round(Tf/dx);
        tbdy = linspace(0,Tf,ntbdy+1)';
        [XX,YY,TT] = meshgrid(xbdy,xbdy,tbdy);
        Ibdy = XX==0 | XX==1 | YY==0 | YY==1 | TT==Tf;
        XX = XX(Ibdy);
        YY = YY(Ibdy);
        TT = TT(Ibdy);
        bdy = [XX(:),YY(:),TT(:)];
         
        xt2  = repel_3d(xt,bdy,@radius_2dTD,rFactor); % Repel by all boundaries
        xt2 = xt2(xt2(:,1)<1 & xt2(:,2)<1 & xt2(:,3)<Tf & xt2(:,1)>0 & xt2(:,2)>0 & xt2(:,3)>0,:);
    case 'lake'
        % generate nodes within the unit circle in the x-y plane
        dx = radius_3d_unif([0,0,0],rFactor);
        R = 1;
        ninit = round(2./dx*10);
        dotmax = round((2./dx)^2*(Tf/dx)*2);
        xt = node_drop_3d([-R R -R R 0 Tf],[ninit,ninit],dotmax,@radius_3d_unif,rFactor);
        % Create boundary nodes
        nbdy = round(2*pi*R/dx);
        dth = 2*pi/nbdy;
        thbdy = -pi:dth:pi-dth;
        xybdy = [R*cos(thbdy)', R*sin(thbdy)'];
        ntbdy = round(Tf/dx);
        tbdy = linspace(0,Tf,ntbdy+1)';
        [XX,TT] = meshgrid(xybdy(:,1),tbdy);
        [YY,~] = meshgrid(xybdy(:,2),tbdy);
        bdy = [XX(:),YY(:),TT(:)];
        % add a bdy at Tf
        xybdyf = node_drop_2d([-R R -R R],ninit,round((2*R/dx)^2*2),@radius_2d,rFactor);
        [thbdyf,rbdyf] = cart2pol(xybdyf(:,1),xybdyf(:,2));
        xybdyf = xybdyf(rbdyf<R-0.5*dx,:);
        tbdyf = Tf*ones(length(xybdyf),1);
        xtbdyf = [xybdyf tbdyf];
        bdy = [bdy; xtbdyf];
        xt2  = repel_3d(xt,bdy,@radius_3d_unif,rFactor); % Repel by all boundaries
        r = sqrt(xt2(:,1).^2+xt2(:,2).^2);
        xt2 = xt2(r<R,:);
    case 'cart'
        % lay down Cartesian grid
        dx = rFactor;
        x = 0:dx:1;
        y = 0:dx:1;
        t = 0:dx:Tf;
        [X,Y,T] = meshgrid(x,y,t);
        X = X(:); Y = Y(:); T = T(:);
        Ibdy = X == 0 | X == 1 | Y == 0 | Y == 1 | T == 0 | T == Tf;
        xt2 = [X(~Ibdy) Y(~Ibdy) T(~Ibdy)];
        bdy = [X(Ibdy) Y(Ibdy) T(Ibdy)];
end



xt = xt2;   % Remove last column with radius information

[ni,~] = size(xt);      % Get total number of interior nodes
[nb,~] = size(bdy);     % Get total number of boundary nodes

xt = [bdy;xt];          % Make xy contain both boundary and interior nodes
% xt = sortrows(xt,3);    % Sort the rows from lowest y (i.e. t) to highest