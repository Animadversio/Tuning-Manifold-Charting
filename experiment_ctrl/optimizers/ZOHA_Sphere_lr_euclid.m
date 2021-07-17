classdef ZOHA_Sphere_lr_euclid < handle
    % Using the Euclidean space mean instead of spherical mean with
    % learning rate
	properties
		dimen   % dimension of input space
        B   % population batch size
        % mu   % scale of the Gaussian distribution to estimate gradient
        lr  % learning rate (step size) of moving along gradient
        sphere_norm 
        tang_codes 
        select_cutoff
        mulist
        mu_init
        mu_final
        
        grad   % estimated gradient
        innerU   % inner random vectors with covariance matrix Id
        outerV   % outer random vectors with covariance matrix H^{-1}, equals innerU @ H^{-1/2}
        xcur   % current base point
        xnew   % new base point

        istep = -1;  % step counter
        counteval = 0
        maximize  % maximize / minimize the function
        rankweight % Switch between using raw score as weight VS use rank weight as score
        rankbasis % Ranking basis or rank weights only
        opts % object to store options for the future need to examine or tune
        
    end

   	methods
   		function self = ZOHA_Sphere_lr_euclid(space_dimen, options)
            self.dimen = space_dimen;  % dimension of input space
            self.parseParameters(options); % parse the options into the initial values of optimizer
            
            self.tang_codes = zeros(self.B, self.dimen);
            self.grad = zeros(1, self.dimen);  % estimated gradient
            self.innerU = zeros(self.B, self.dimen);  % inner random vectors with covariance matrix Id
            self.outerV = zeros(self.B, self.dimen);  % outer random vectors with covariance matrix H^{-1}, equals innerU @ H^{-1/2}
            self.xcur = zeros(1, self.dimen); % current base point
            self.xnew = zeros(1, self.dimen); % new base point
            
   		end % of initialization
        
        function parseParameters(self, opts)
        	if ~isfield(opts, "population_size"), opts.population_size = 40; end
    		if ~isfield(opts, "select_cutoff"), opts.select_cutoff = opts.population_size / 2; end
    		%if ~isfield(opts, "mu"), opts.mu = 0.005; end
    		if ~isfield(opts, "lr"), opts.lr = 2; end
    		if ~isfield(opts, "sphere_norm"), opts.sphere_norm = 300; end
    		if ~isfield(opts, "maximize"), opts.maximize = true; end
    		if ~isfield(opts, "rankweight"), opts.rankweight = true; end
            if ~isfield(opts, "rankbasis"), opts.rankbasis = false; end
            if ~isfield(opts, "mu_init"), opts.mu_init = 0.01; end  % need to check
            if ~isfield(opts, "mu_final"), opts.mu_final = 0.002; end
            if ~isfield(opts, "indegree"), opts.indegree = false; end
            self.B = opts.population_size;  % population batch size
            self.select_cutoff = floor(opts.select_cutoff);
            %self.mu = opts.mu;  % scale of the Gaussian distribution to estimate gradient
            if ~ opts.indegree % mu in absolute units 
                self.mu_init = opts.mu_init; 
                self.mu_final = opts.mu_final; 
            else % mu in degrees
                self.mu_init = opts.mu_init / 180 * pi / sqrt(self.dimen); 
                self.mu_final = opts.mu_final / 180 * pi / sqrt(self.dimen); 
            end
            self.lr = opts.lr;  % learning rate (step size) of moving along gradient
            self.sphere_norm = opts.sphere_norm;
            self.maximize = opts.maximize;  % maximize / minimize the function
            self.rankweight = opts.rankweight; % Switch between using raw score as weight VS use rank weight as score
            self.rankbasis = opts.rankbasis; % whether include basis in the ranking comparison.
            self.mulist = []; % do initialize use the functio lr_schedule
            fprintf("\nSpherical Space dimension: %d, Population size: %d, Optimization Parameters:\n Exploration: %.3f-%.3f\n Learning rate: %.3f\n Sphere Radius: %.3f\n",...
               self.dimen, self.B, self.mu_init, self.mu_final, self.lr, self.sphere_norm)
            if self.rankweight
            	fprintf("Using rank weight, selection size: %d\n", self.select_cutoff)
            end
            self.opts = opts; % save a copy of opts with default value updated. Easy for printing.
            % Parameter Checking 
            ExpectExplAng1 = (sqrt(self.dimen) * self.mu_init) / pi * 180; % Expected angular distance between sample and basis 
            ExpectExplAng2 = (sqrt(self.dimen) * self.mu_final) / pi * 180; % Expected angular distance between sample and basis 
            fprintf("Expected angular exploration length %.1f - %.1f deg\n", ExpectExplAng1, ExpectExplAng2)
            if ExpectExplAng1 > 90 || ExpectExplAng2 > 90
                warning("Estimated exploration range too large! Destined to fail! Check parameters!\n")
            end
            if self.rankweight
                weights = rankweight(self.B, self.select_cutoff);
                ExpectStepSize1 = sqrt(self.dimen * sum(weights.^2)) * self.mu_init * self.lr / pi * 180;
                ExpectStepSize2 = sqrt(self.dimen * sum(weights.^2)) * self.mu_final * self.lr / pi * 180;
                % Expected angular distance of gradient step.  
                fprintf("Estimated angular step size %.1f - %.1f deg\n", ExpectStepSize1, ExpectStepSize2)
                if ExpectStepSize1 > 90 || ExpectStepSize2 > 90
                    warning("Estimated step size too large! Destined to fail! Check parameters!\n")
	            end
            end
        end
        
        function lr_schedule(self,gen_total,mode)
            if nargin == 2
                mode = "exp";
            end
            switch mode
                case "lin"
                self.mulist = linspace(self.mu_init, self.mu_final, gen_total);
                case "exp"
                self.mulist = logspace(log10(self.mu_init), log10(self.mu_final), gen_total);
                case "inv"
                self.mulist = 15 + 1./(0.0017 * [1:gen_total] + 0.0146);
                self.opts.mu_init = self.mulist(1);
                self.opts.mu_final = self.mulist(end);
                self.mulist = self.mulist / 180 * pi / sqrt(self.dimen);
                self.mu_init = self.mulist(1); self.mu_final = self.mulist(end);
            end
        end
        
   		function [new_samples, new_ids] =  doScoring(self, codes, scores, maximize, TrialRecord)
   		N = self.dimen;
        fprintf('Gen %d max score %.3f, mean %.3f, std %.3f\n',...
                self.istep, max(scores),mean(scores),std(scores) )
        if self.istep == -1
            % Population Initialization: if without initialization, the first xmean is evaluated from weighted average all the natural images
            fprintf('First generation\n')
            self.xcur = codes(1, :);
            if self.rankweight == false % use the score difference as weight
                % B normalizer should go here larger cohort of codes gives more estimates 
                weights = (scores - scores(1)) / self.B; % / self.mu 
            else  % use a function of rank as weight, not really gradient. 
                if self.maximize == false % note for weighted recombination, the maximization flag is here. 
                    [~,code_rank]=ismember(scores, sort(scores,'ascend')); % find rank of ascending order
                else
                    [~,code_rank]=ismember(scores, sort(scores,'descend')); % find rank of descending order 
                end
                % Note the weights here are internally normalized s.t. sum up to 1, no need to normalize more. 
                raw_weights = rankweight(length(code_rank)); 
                weights = raw_weights(code_rank); % map the rank to the corresponding weight of recombination
                % Consider the basis in our rank! but the weight will be wasted as we don't use it. 
            end
            w_mean = weights * codes; % mean in the euclidean space
            self.xnew = w_mean / norm(w_mean) * self.sphere_norm; % project it back to shell. 
        else
            % self.xcur = self.xnew % should be same as following line
            self.xcur = codes(1, :);
            if self.rankweight == false % use the score difference as weight
                % B normalizer should go here larger cohort of codes gives more estimates 
                weights = (scores(1:end) - scores(1)) / self.B; % / self.mu 
            else  % use a function of rank as weight, not really gradient. 
                if ~ self.rankbasis % if false, then exclude the first basis vector from rank (thus it receive no weights.)
                    rankedscore = scores(2:end);
                else
                    rankedscore = scores;
                end
                if self.maximize == false % note for weighted recombination, the maximization flag is here. 
                    [~,code_rank]=ismember(rankedscore, sort(rankedscore,'ascend')); % find rank of ascending order
                else
                    [~,code_rank]=ismember(rankedscore, sort(rankedscore,'descend')); % find rank of descending order 
                end
                % Note the weights here are internally normalized s.t. sum up to 1, no need to normalize more. 
                raw_weights = rankweight(length(code_rank)); 
                weights = raw_weights(code_rank); % map the rank to the corresponding weight of recombination
                % Consider the basis in our rank! but the weight will be wasted as we don't use it. 
                if ~ self.rankbasis
                    weights = [0, weights]; % the weight of the basis vector will do nothing! as the deviation will be nothing
                end
            end
            % estimate gradient from the codes and scores
            % assume self.weights is a row vector 
            w_mean = weights * codes; % mean in the euclidean space
            w_mean = w_mean / norm(w_mean) * self.sphere_norm; % project it back to shell. 
            self.xnew = SLERP(self.xcur, w_mean, self.lr); % use lr to spherical extrapolate
            fprintf("Step size %.3f, multip learning rate %.3f, ",  ang_dist(self.xcur, self.xnew), ang_dist(self.xcur, self.xnew) * self.lr);
            ang_basis_to_samp = ang_dist(codes, self.xnew); 
            fprintf("New basis ang to last samples mean %.3f(%.3f), min %.3f\n",  mean(ang_basis_to_samp), std(ang_basis_to_samp), min(ang_basis_to_samp));
%             HAgrad = weights * self.tang_codes;
%             fprintf("Estimated Gradient Norm %f\n",  norm(HAgrad))
%             % determine the moving direction and move. 
%             if (~self.maximize) && (~self.rankweight), mov_sign = -1; else, mov_sign = 1; end
%             self.xnew = ExpMap(self.xcur, mov_sign * self.lr * HAgrad);
        end
        % Generate new sample by sampling from Gaussian distribution
        self.tang_codes = zeros(self.B, N);  % Tangent vectors of exploration
        new_samples = zeros(self.B + 1, N);
        self.innerU = randn(self.B, N);  % Isotropic gaussian distributions
        self.outerV = self.innerU; % H^{-1/2}U, more transform could be applied here! 
        self.outerV = self.outerV - (self.outerV * self.xnew') * self.xnew / norm(self.xnew)^2; % orthogonal projection to xnew's tangent plane.
        
        new_samples(1, :) = self.xnew;
        self.tang_codes = self.mulist(self.istep + 2) * self.outerV; % m + sig * Normal(0,C)
        new_samples(2:end, :) = ExpMap(self.xnew, self.tang_codes); 
        fprintf("Current Exploration %.1f deg\n", self.mulist(self.istep + 2) * sqrt(self.dimen - 1) / pi * 180)
        new_ids = [];
        for k = 1:size(new_samples,1)
            new_ids = [new_ids, sprintf("gen%03d_%06d", self.istep+1, self.counteval)]; 
            self.counteval = self.counteval + 1;
        end
        self.istep = self.istep + 1;
        new_samples = renormalize(new_samples, self.sphere_norm);
        end
    end
end