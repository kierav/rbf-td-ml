function W = RBF_TD_PHS_AdvDiff2D(x,y,t,xc,yc,tc,alpha,beta,tscale,m,d)

% Input parameters
%   x,y,t     Column vectors with the x-,y- and t-coordinates within the stencil
%   xc,yc,tc   Stencil center locations
%   alpha   Local wave speed
%   m       Order of PHS
%   d       Polynomial order to attach
%   tscale  Scaling of the time dimension
% Output parameters
%     W     RBF-FD weights for u_t+alpha_x*u_x + alpha_y*u_y = beta*(u_xx
%     + u_yy)

if length(beta) < 2
    beta = [beta beta];
end

%Shift nodes
xn = x(1); yn = y(1); tn = t(1);
x = x - xn; 
y = y - yn;
t = t - tn;
xc = xc - xn;
yc = yc - yn;
tc = tc - tn;
n = length(x);
D1 = sqrt((x-x').^2+(y-y').^2+(tscale*(t-t')).^2);

[x1,x2] = meshgrid(xc,x);           % Calculate the weights for
[y1,y2] = meshgrid(yc,y);
[t1,t2] = meshgrid(tc,t);           % the p tentative evaluation 
xd   = x1-x2;  yd = y1-y2;   td = t1-t2;       % locations
D2 = sqrt(xd.^2 + yd.^2 + (tscale*td).^2);

fi  =@(r) r.^m;     %PHS; Radial function
Lfi =@(r,x,y,t) -m*(tscale^2*t+alpha(1)*x+alpha(2)*y-beta(1)-beta(2)).*r.^(m-2)...
    + m*(m-2)*(beta(1)*x.^2+beta(2)*y.^2).*r.^(m-4) ; %Operator acting on PHS

A11 = fi(D1);                % RBF matrix within stencil
L1 = Lfi(D2,xd,yd,td); 

% appended polynomials
% P = d; % degree of polynomial terms
% x_exp = []; y_exp; t_exp = [];
% for deg = 0:P
%     x_exp = [x_exp,deg:-1:0];
%     t_exp = [t_exp,0:1:deg];
% end 
% nP = length(x_exp);

if d == -1                     % Special case; no polynomial terms,
    A = A11;  L = L1;    
else                        % Include polynomials and matching constraints
    X   =  x(:,ones(1,d+1));  X(:,1) = 1;  X = cumprod(X,2);
    Y   =  y(:,ones(1,d+1));  Y(:,1) = 1;  Y = cumprod(Y,2);
    Z   =  t(:,ones(1,d+1));  Z(:,1) = 1;  Z = cumprod(Z,2);
    np  = (d+3)*(d+2)*(d+1)/6;        % Number of polynomial terms
    XYZ = zeros(n,np); col = 0;       % Assemble polynomial matrix block
    
    nc = length(xc);
    L2 = zeros(np,nc);
    
    for kt = 0:d
        for kz = 0:kt
            for ky = 0:kt-kz
                col = col+1;
                kx = kt-ky-kz;
                XYZ(:,col) = X(:,kx+1).*Y(:,ky+1).*Z(:,kz+1);
                L2(col,:) = -(kz*xc'.^kx.*yc'.^ky.*tc'.^(max(0,kz-1)) + ...
                            alpha(1)*kx*xc'.^(max(0,kx-1)).*yc'.^ky.*tc'.^kz + ...
                            alpha(2)*ky*xc'.^kx.*yc'.^(max(0,ky-1)).*tc'.^kz) + ...
                            beta(1)*kx*(kx-1)*xc'.^(max(0,kx-2)).*yc'.^ky.*tc'.^kz + ...
                            beta(2)*ky*(ky-1)*xc'.^kx.*yc'.^(max(0,ky-2)).*tc'.^kz;
             end
        end
    end

%     fx = (x_exp(ones(1,nc),:).*xc(:,ones(1,nP)).^(max(0,x_exp(ones(1,nc),:)-1)).*tc(:,ones(1,nP)).^t_exp(ones(1,nc),:))';
%     ft = (t_exp(ones(1,nc),:).*tc(:,ones(1,nP)).^(max(0,t_exp(ones(1,nc),:)-1)).*xc(:,ones(1,nP)).^x_exp(ones(1,nc),:))';
%     L2 = -(ft+alpha*fx);
    
%     if(d>=1) L2(2)=alpha; L2(3)=1; end
    
    L = [L1;L2];     % Assemble the linear system to be solved
    
    A = [A11,XYZ;XYZ',zeros(col,col)];
    
end

w = A\L;
W = w(1:n,:);