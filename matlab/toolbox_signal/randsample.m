function y = randsample(n, k, replace, w)
%RANDSAMPLE Random sample, with or without replacement.
%   Y = RANDSAMPLE(N,K) returns Y as a vector of K values sampled uniformly
%   at random, without replacement, from the integers 1:N.
%
%   Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at
%   random, without replacement, from the values in the vector POPULATION.
%
%   Y = RANDSAMPLE(...,REPLACE) returns a sample taken with replacement if
%   REPLACE is true, or without replacement if REPLACE is false (the default).
%
%   Y = RANDSAMPLE(...,true,W) returns a weighted sample, using positive
%   weights W, taken with replacement.  W is often a vector of probabilities.
%   This function does not support weighted sampling without replacement.
%
%   Example:  Generate a random sequence of the characters ACGT, with
%   replacement, according to specified probabilities.
%
%      R = randsample('ACGT',48,true,[0.15 0.35 0.35 0.15])
%
%   See also RAND, RANDPERM.

%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 1.1.4.3 $  $Date: 2008/12/01 08:09:34 $

if nargin < 2
    error('stats:randsample:TooFewInputs','Requires two input arguments.');
elseif numel(n) == 1
    population = [];
else
    population = n;
    n = numel(population);
    if length(population)~=n
       error('stats:randsample:BadPopulation','POPULATION must be a vector.');
    end
end

if nargin < 3
    replace = false;
end

if nargin < 4
    w = [];
elseif ~isempty(w)
    if length(w) ~= n
        if isempty(population)
            error('stats:randsample:InputSizeMismatch',...
                  'W must have length equal to N.');
        else
            error('stats:randsample:InputSizeMismatch',...
                  'W must have the same length as the population.');
        end
    else
        p = w(:)' / sum(w);
    end
end

switch replace

% Sample with replacement
case {true, 'true', 1}
    if isempty(w)
        y = ceil(n .* rand(k,1));
    else
        [dum, y] = histc(rand(k,1),[0 cumsum(p)]);
    end

% Sample without replacement
case {false, 'false', 0}
    if k > n
        if isempty(population)
            error('stats:randsample:SampleTooLarge',...
        'K must be less than or equal to N for sampling without replacement.');
        else
            error('stats:randsample:SampleTooLarge',...
                  'K must be less than or equal to the population size.');
        end
    end

    if isempty(w)
        % If the sample is a sizeable fraction of the population,
        % just randomize the whole population (which involves a full
        % sort of n random values), and take the first k.
        if 4*k > n
            rp = randperm(n);
            y = rp(1:k);

        % If the sample is a small fraction of the population, a full sort
        % is wasteful.  Repeatedly sample with replacement until there are
        % k unique values.
        else
            x = zeros(1,n); % flags
            sumx = 0;
            while sumx < k
                x(ceil(n * rand(1,k-sumx))) = 1; % sample w/replacement
                sumx = sum(x); % count how many unique elements so far
            end
            y = find(x > 0);
            y = y(randperm(k));
        end
    else
        error('stats:randsample:NoWeighting',...
              'Weighted sampling without replacement is not supported.');
    end
otherwise
    error('stats:randsample:BadReplaceValue',...
          'REPLACE must be either true or false.');
end

if ~isempty(population)
    y = population(y);
else
    y = y(:);
end
