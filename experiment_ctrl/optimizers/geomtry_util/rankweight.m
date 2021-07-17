function weights = rankweight(lambda_, mu)
    % Rank weight adapted from CMA-ES code
    % mu is the cut off number, how many samples will be kept while `lambda_ - mu` will be ignore 
    if nargin == 1
        mu = lambda_ / 2;  % number of parents/points for recombination
        %  Defaultly Select half the population size as parents
    end
    weights = zeros(1, int32(lambda_));
    mu_int = floor(mu);
    weights(1:mu_int) = log(mu + 1 / 2) - log([1 : mu_int]);  % muXone array for weighted recombination
    weights = weights / sum(weights);
return 