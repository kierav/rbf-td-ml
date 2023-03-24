function [X,y,xtU,u,errL2Tf,errLinfTf,time] = RBF_TD_Wave2D_func(alpha,bc,u_sol,Tf,rFactor,nodetype,useMLweights,useL2,plotFigs,modelname)
% RBF-TD solver for the 2-D wave problem 
%   utt = alpha*(u_xx + u_yy)
%
% Input parameters:
%   alpha       wavespeed 
%   bc          boundary condition (p - periodic)
%   u_sol       function in (x,y,t) for analytical solution 
%   Tf          final time
%   rFactor     scaling factor for nodeset
%   nodetype    type of nodeset to generate (quasiunif, cart, lake or
%                   vardensity)
%   useMLweights    1 if use ML model for weight optimization, else 0
%   useL2       1 if use L2 for weight optimization, else 0
%   plotfigs    1 if plot solution figures
%   modelname   filename for ML model
%
% Output parameters:
%   X           normalized candidate weight sets (input to weight optimization)
%   y           combination weights (output of weight optimization)
%   xtU         (x,y,t) nodeset
%   u           numerical solution on nodes specified by xtU
%   errL2Tf     L2 error along bdyf
%   errLinfTf   Linf error along bdyf
%   time        clock time to perform weight optimization

    % Parameters for RBF stencils
    n = 70;               % Size of weight stencils
    p = 9;                % Number of candidate stencils to generate
    constp = 5;           % Subset of candidate stencils to preselect 
    d = 4;                % Degree of poynomials to attach
    m = 3;                % Degree of PHS

    % Generate node set
    xt = create_nodes_2d(rFactor,nodetype,Tf);

    % Reorder node set
    if bc == 'p'
        xtU = xt(xt(:,2)<1 & xt(:,1)<1 & xt(:,3)<=Tf,:);
    else
        xtU = xt(xt(:,2)<=1 & xt(:,1)<=1 & xt(:,3)<=Tf,:);
    end
    
    xtU = sortrows(xtU,3);
    bdyf = abs(xtU(:,3)-Tf)<1e-10;
    [nt,~] = size(xtU);    

    % define boundary
    if strcmp(nodetype,'quasiunif') && bc ~= 'p'
        bdy = xtU(:,1) < 1e-8 | xtU(:,1) > 1-1e-8 | xtU(:,2) < 1e-8 | xtU(:,2) > 1-1e-8;
    elseif strcmp(nodetype,'lake')
        r = sqrt(xtU(:,1).^2+xtU(:,2).^2);
        bdy = r>1-1e-8;
    else
        bdy = zeros(nt,1);
    end

    % Find nearest neighbors
    if bc == 'p'
        xtUtiled = [xtU; xtU(:,1)-1 xtU(:,2)+1 xtU(:,3); xtU(:,1) xtU(:,2)+1 xtU(:,3);
            xtU(:,1)+1 xtU(:,2)+1 xtU(:,3); xtU(:,1)-1 xtU(:,2) xtU(:,3);
            xtU(:,1)+1 xtU(:,2) xtU(:,3); xtU(:,1)-1 xtU(:,2)-1 xtU(:,3);
            xtU(:,1) xtU(:,2)-1 xtU(:,3); xtU(:,1)+1 xtU(:,2)-1 xtU(:,3);];
        splits = 6;
        lengthIDX = ceil(nt/splits);
        buffer = 2*round((nt).^(1/3)*n);
        IDX = [];
        for i = 1:splits
            indL = max(1,(i-1)*lengthIDX-buffer);
            indR = min(nt,i*lengthIDX+buffer);
            [Idx,~] = knnsearch(xtUtiled,xtU((i-1)*lengthIDX+1:min(i*lengthIDX,nt),:),'K',12*n);
            IDX = [IDX;Idx];
        end
        IDX = mod(IDX,nt);
        IDX(IDX == 0) = nt;
    else
        [IDX,~] = knnsearch(xtU,xtU,'K',20*n);
    end


    u       = zeros(nt,1);  % Vector to hold computed solution
    weights = zeros(nt,n);  % Array to hold computed weights for all stencils
    W       = zeros(n,constp,nt);% Array to hold all the temporary weight sets
            % Dimensions run over: (1) Weights, (2) Eval. points, (3) Stencils
    ind2    = zeros(nt,n);  % Array to hold lower neighbors to node k
    iBC     = zeros(nt,1);  % Array to mark nodes with given boundary values
    Wscaled = zeros(n,constp,nt);   % Array to hold normalized weight sets

    for k = 1:nt            % Sequential loop over all nodes
                            % Create the weight sets to be optimized in a later
                            % parfor loop
        if xtU(k,3) < 3*rFactor || (bdy(k) ==1 && bc~='d') % Find k-values for bdy nodes
            u(k)   = u_sol(xtU(k,1),xtU(k,2),xtU(k,3));     % Assign given boundary value
            iBC(k) = true;                      % Mark that this has been done
        elseif bdy(k) == 1 && bc=='d'
            u(k) = 0;
            iBC(k) = true;
        else
            % Locate the n nearest neighbors below node k
            if alpha>1
                sigma = 1/alpha;
                scaling = 1/sigma;
                %Finding neighbors
                [ind,~]=knnsearch(xtU(xtU(:,3)<xtU(k,3),:),xtU(k,:),'K',n-1,'Distance','mahalanobis','Cov',[1 0 0; 0 1 0; 0 0 sigma^2]);
                ind = [k,ind];
            else
                scaling = 1;
                ind = IDX(k,:);
                ind = [ind(1),ind(xtU(ind,3)<xtU(k,3))];
                if length(ind)>=n
                    ind = ind(1:n);
                else
                    fprintf('Not enough neighbours\n');
                end
            end

            ind2(k,:) = ind;

            x = xtU(ind,1);  y = xtU(ind,2);    t = xtU(ind,3);  % Get stencil node locations x,t

            if bc == 'p'
                % deal with periodic boundaries
                xl = 0; xr = 1;
                x = x - x(1); y = y - y(1);
                x(x>(xr-xl)/2) = x(x>(xr-xl)/2)-(xr-xl);
                x(x<-(xr-xl)/2) = x(x<-(xr-xl)/2)+(xr+xl);
                y(y>(xr-xl)/2) = y(y>(xr-xl)/2)-(xr-xl);
                y(y<-(xr-xl)/2) = y(y<-(xr-xl)/2)+(xr+xl);
                x = x + x(1); y = y + y(1);
            end
            
            % Center points
            xe      = x(1)+0.5*(x(2:p+1)-x(1));
            ye      = y(1)+0.5*(y(2:p+1)-y(1));
            te      = t(1)+0.5*(t(2:p+1)-t(1));     

            % Calculate all the weights
            w = RBF_TD_PHS_Wave2D(x,y,t,xe,ye,te,alpha,scaling,m,d);
            
            % Select subset of "best" candidate stencils
            if p > constp
                wmax = max(abs(w),[],1);
                [~,I] = sort(wmax-abs(w(1,:)));
                I = I(1:constp);
                I = sort(I);
                w = w(:,I);
            end
            W(:,:,k) = w;

            % Normalize stencils
            fact = 1./w(1,:);
            W2 = bsxfun(@times,w,fact);
            W2(:,2:constp) = bsxfun(@minus,W2(:,2:constp),W2(:,1));
            Wscaled(:,:,k) = W2;
        end
    end
  
    X = Wscaled(:,:,~iBC);

    % Optimize weights over candidate stencils
    [dvs,~,time] = optimizeWeights(useMLweights,useL2,Wscaled,iBC,modelname);
    
    y = dvs(:,~iBC);

    % Apply all the stencils for the actual time stepping
    for k = 1:nt
        if ~iBC(k)
            dv = dvs(:,k);
            wf = Wscaled(:,1,k)-Wscaled(:,2:constp,k)*dv;
            weights(k,:) = wf';     % Store weights for actual time stepping   
            u(k) = -weights(k,2:n)*u(ind2(k,2:n))/weights(k,1); % Apply stencil
        end
    end
    
    %Error in solution
    er = u - u_sol(xtU(:,1),xtU(:,2),xtU(:,3));
    errL2 = norm(er,2)/sqrt(length(xtU(:,1)));
    errLinf = norm(er,inf);
    
    %Error at final time
    errTf = u(bdyf)-u_sol(xtU(bdyf,1),xtU(bdyf,2),xtU(bdyf,3));
    errL2Tf = norm(errTf,2)/sqrt(length(xtU(bdyf,1)));
    errLinfTf = norm(errTf,inf);

    if plotFigs
        % interpolate data to a grid for visualization
        Fnum = scatteredInterpolant(xtU,real(u));
        FnumTf = scatteredInterpolant(xtU(bdyf,1:2),real(u(bdyf)));
        
        % create visualization grid
        xx = linspace(0,1,25);
        [Xq,Yq,Tq] = meshgrid(xx,xx,xx);
        Uq = Fnum(Xq,Yq,Tq);
        Utrueq = real(u_sol(Xq,Yq,Tq));
        
        % plot volumetric slices
        xslice = [0.1,0.5,1];
        yslice = [];
        zslice = 0.5;
        figure(1)
        slice(Xq,Yq,Tq,Uq,xslice,yslice,zslice);
        title('Numeric solution through time')
        caxis([-1 1])
        figure(2)
        slice(Xq,Yq,Tq,Utrueq-Uq,xslice,yslice,zslice);
        title('Error in numeric solution through time')
        colorbar
        
        % visualize final solution
        xx = linspace(min(xtU(:,1)),max(xtU(:,1)),25);
        [Xg,Yg] = meshgrid(xx,xx);
        Ug = FnumTf(Xg,Yg);
        if strcmp(nodetype,'lake') 
            r = sqrt(Xg.^2+Yg.^2);
            Ug(r>1) = NaN;
        end
        figure(3)
        mesh(Xg,Yg,Ug);
        title('Numeric solution at Tf')
        figure(4)
        mesh(Xg,Yg,Ug-real(u_sol(Xg,Yg,Tf)))
        title('Error in numeric solution at Tf');
    end
end