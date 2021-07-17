% Initially code in 2019. 
% Adapted from the python version of CMAES optimizer by Binxu.
% 
% Binxu & Carlos.
% Dec. 2020
classdef CMAES_simple < handle
    
    properties
        codes
        scores
        sigma
        N
        istep
        mu
        mueff
        lambda
        damps
        chiN
        weights
        c1
        cs
        cc
        ps
        pc
        A
        Ainv
        Aupdate_freq = 10 ;
        randz
        counteval
        eigeneval
        update_crit
        xmean
        init_x
        opts
    end % of properties
    
    
    methods
        
        function obj = CMAES_simple(space_dimen, init_x, options) 
            % 2nd should be [] if we do not use init_x parameter. 
            % if we want to tune this algorithm, we may need to send in a
            % structure containing some initial parameters? like `sigma`
            % object instantiation and parameter initialization 
            % obj.codes = codes; % not used..? Actually the codes are set in the first run
            obj.N = space_dimen; % size(codes,2);
            if ~isfield(options, "popsize")
            obj.lambda = 4 + floor(3 * log2(obj.N));  % population size, offspring number
            % the relation between dimension and population size.
            else
            obj.lambda = options.popsize;
            end
            obj.randz = randn(obj.lambda, obj.N);  % Not Used.... initialize at the end of 1st generation
            if ~isfield(options, "mu")
            obj.mu = obj.lambda / 2;  % number of parents/points for recombination
            else
            obj.mu = options.mu;
            end
            %  Select half the population size as parents
            obj.weights = log(obj.mu + 1/2) - log( 1:floor(obj.mu));  % muXone array for weighted recombination
            obj.mu = floor(obj.mu);
            obj.weights = obj.weights / sum(obj.weights);  % normalize recombination weights array
            obj.mueff = sum(obj.weights).^ 2 / sum(obj.weights.^ 2);
            
            obj.cc = 4 / (obj.N + 4);
            obj.cs = sqrt(obj.mueff) / (sqrt(obj.mueff) + sqrt(obj.N));
            obj.c1 = 2 / (obj.N + sqrt(2))^2;
            obj.damps = 1 + obj.cs + 2 * max(0, sqrt((obj.mueff - 1) / (obj.N + 1)) - 1);  % damping for sigma usually close to 1
            obj.chiN = sqrt(obj.N) * (1 - 1 / (4 * obj.N) + 1 / (21 * obj.N^2));
            % expectation of ||N(0,I)|| == norm(randn(N,1)) in 1/N expansion formula
            
            obj.pc = zeros(1, obj.N);
            obj.ps = zeros(1, obj.N);
            obj.A = eye(obj.N, obj.N);
            obj.Ainv = eye(obj.N, obj.N);
            obj.eigeneval=0;
            obj.counteval=0;
            
            
            % if init_x is set in TrialRecord, use it, it will become the first
            if ~isempty(init_x)
                obj.init_x = init_x;
            else
                obj.init_x = [];
            end
            % xmean in 2nd generation
            obj.xmean = zeros(1, obj.N); % Not used.
            if ~isfield(options, "init_sigma"), options.init_sigma = 3;end
            if ~isfield(options, "Aupdate_freq"), options.Aupdate_freq = 10; end
            obj.Aupdate_freq = options.Aupdate_freq;
            obj.update_crit = obj.lambda * obj.Aupdate_freq;% / obj.c1 / obj.N ;
            obj.sigma = options.init_sigma; 
            obj.istep = -1;
            obj.opts = options;
        end % of initialization
        
        
        function [new_samples, new_ids, TrialRecord] =  doScoring(obj,codes,scores,maximize,TrialRecord)
            
            obj.codes = codes;
            
            % Sort by fitness and compute weighted mean into xmean
            if ~maximize
                [sorted_score, code_sort_index] = sort(scores);  % add - operator it will do maximization.
            else
                [sorted_score, code_sort_index] = sort(scores, 'descend');
            end
            % TODO: maybe print the sorted score?
            % disp(sorted_score')
%             fprintf("max score %.3f, mean %.3f (std %.3f)",...
%                 max(sorted_score),mean(sorted_score),std(sorted_score) )
            
            if obj.istep == -1 % if first step
%                 fprintf('is first gen\n');
                % Population Initialization: if without initialization, the first obj.xmean is evaluated from weighted average all the natural images
                if isempty(obj.init_x) ||  size(obj.init_x, 1)~=1
                    if obj.mu<=length(scores)
                        obj.xmean = obj.weights * obj.codes(code_sort_index(1:obj.mu), :);
                    else % if ever the obj.mu (selected population size) larger than the initial population size (init_population is 1 ?)
                        tmpmu = max(floor(length(scores)/2), 1); % make sure at least one element!
                        obj.xmean = obj.weights(1:tmpmu) * obj.codes(code_sort_index(1:tmpmu), :) / sum(obj.weights(1:tmpmu));
                    end
                else
                    obj.xmean = obj.init_x;
                end
                
                
            else % if not first step
                
%                 fprintf('not first gen\n');
%                 xold = obj.xmean;
                % Weighted recombination, move the mean value
                obj.xmean = obj.weights * obj.codes(code_sort_index(1:obj.mu), :);
                
                % Cumulation: Update evolution paths
                randzw = obj.weights * obj.randz(code_sort_index(1:obj.mu), :); % Note, use the obj.randz saved from last generation
                obj.ps = (1 - obj.cs) * obj.ps + sqrt(obj.cs * (2 - obj.cs) * obj.mueff) * randzw;
                obj.pc = (1 - obj.cc) * obj.pc + sqrt(obj.cc * (2 - obj.cc) * obj.mueff) * randzw * obj.A;
                
                % Adapt step size sigma
                obj.sigma =  obj.sigma * exp((obj.cs / obj.damps) * (norm(obj.ps) / obj.chiN - 1));
                
                %     if std(scores) < 5 && obj.istep > 50
                %         sigma = sigma + sigma.*0.10 ;
                %     end
                %     if obj.istep > 5
                % %     if all(~diff(TrialRecord.User.respOverGen(end-4:end))), sigma = 3 ; end
                %     end
                
                fprintf("Step %d, sigma: %0.2e, Scores\n",obj.istep, obj.sigma)
                
                % Update the factors of C,  (diagonalization)
                if obj.counteval - obj.eigeneval > obj.update_crit  % to achieve O(N ^ 2)
                    obj.eigeneval = obj.counteval;
                    tic;
                    v = obj.pc * obj.Ainv;
                    normv = v * v';
                    obj.A = sqrt(1-obj.c1) * obj.A + sqrt(1-obj.c1)/normv*(sqrt(1+normv*obj.c1/(1-obj.c1))-1) * v * obj.pc';
                    obj.Ainv = 1/sqrt(1-obj.c1) * obj.Ainv - 1/sqrt(1-obj.c1)/normv*(1-1/sqrt(1+normv*obj.c1/(1-obj.c1))) * obj.Ainv* (v'*v);
                    fprintf("obj.A, obj.Ainv update! Time cost: %.2f s\n", toc)
                end
                
            end % of first step
            
            % Generate new sample by sampling from Multivar Gaussian distribution
            new_samples = zeros(obj.lambda, obj.N);
            new_ids = [];
            obj.randz = randn(obj.lambda, obj.N);  % save the random number for generating the code.
            % For optimization path update in the next generation.
            
            for k = 1:obj.lambda
                new_samples(k, :) = obj.xmean + obj.sigma * (obj.randz(k, :) * obj.A);  % m + sig * Normal(0,C)
                % Clever way to generate multivariate gaussian!!
                % Stretch the guassian hyperspher with D and transform the
                % ellipsoid by B mat linear transform between coordinates
                new_ids = [new_ids, sprintf("gen%03d_%06d",obj.istep+1, obj.counteval)];
               
                % FIXME obj.A little inconsistent with the naming at line 173/175/305/307 esp. for gen000 code
                % assign id to newly generated images. These will be used as file names at 2nd round
                obj.counteval = obj.counteval + 1;
            end
            new_ids = cellstr(new_ids); % Important to save the cell format string array!
            
            % self.sigma, self.obj.A, self.obj.Ainv, self.obj.obj.ps, self.obj.pc = sigma, obj.A, obj.Ainv, obj.obj.ps, obj.pc,
            obj.istep = obj.istep + 1;
%             fprintf('step is now %d\n',obj.istep);
            TrialRecord.User.sigma = obj.sigma ; 
            

        end % of do scoring

    end % of properties
    
end % of classdef
