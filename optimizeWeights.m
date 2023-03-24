function [dvs,dvs_diff,time] = optimizeWeights(useMLWeights,useL2,Wscaled,iBC,modelname)
% optimizeWeights generates the stencil combination weights for RBF-TD
% 
% --- Input parameters ---
%   useMLWeights - 1 or 0 
%   useL2        - 1 or 0
%   Wscaled      - scaled stencil weights
%   iBC          - indices that do not have weights associated (IC or BC nodes)
%   modelname    - string to name the ML model 
% --- Output parameters ---
%   dvs          - stencil combination weights
%   dvs_diff     - when using the ML L2-L1 model, these are the ML updates
%                  which are used to adjust the L2-min generated stencil weights
%   time         - clock time to generate the optimized weights
    
    n = size(Wscaled,1); constP = size(Wscaled,2); nt = size(Wscaled,3);
    dvs = zeros(constP-1,nt);
    dvs_diff = zeros(constP-1,nt);
    
    if useMLWeights     % generate stencil combination weights using a trained ML model
        X = Wscaled(:,:,~iBC);
        save('test_data.mat','X','-v7.3');  % save data for input to ML model
        save('model_data.mat','modelname'); % model name for saving
        system('python model2.py');         % run model on saved data
%         source /curc/sw/anaconda/default; conda activate learn; 
        load('predicted_data.mat','y_pred','time'); % load model output
        dvs(:,~iBC) = y_pred';
        if useL2    % add model output to L2-minimization result to get final stencil combination
            dvs_diff = dvs;
            tic
            for k = 1:nt
                if ~iBC(k)
                    dvl2 = lsqminnorm(Wscaled(2:n,2:constP,k),Wscaled(2:n,1,k));
                    dvs(:,k) = dvl2 + dvs(:,k);
                end
            end
            time = time + toc;
        end
    else
        if useL2    % generate stencil combination weights using L2-minimization
            tic
            for k = 1:nt       
                if ~iBC(k)          % If node not a boundary node, then optimize
                    dv = lsqminnorm(Wscaled(2:n,2:constP,k),Wscaled(2:n,1,k));
                    dvs(:,k) = dv;    
                end
            end  
            time = toc;
        else        % generate stencil combination weights using L1-minimization
            tic
            for k = 1:nt        
                if ~iBC(k)          % If node not a boundary node, then optimize
                    dv = l1decode_pd(zeros(constP-1,1)/(constP-1),Wscaled(2:n,2:constP,k),[],Wscaled(2:n,1,k));
                    dvs(:,k) = dv;    
                end
            end  
            time = toc;
        end
    end
        
end

