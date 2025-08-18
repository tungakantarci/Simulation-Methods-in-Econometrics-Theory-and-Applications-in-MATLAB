function [H,Z] = halton(N, dimensions, draws, varargin)
%HALTON Generate Halton sequences and transform to standard normal.
% Version 1.0, Updated August 2025
% Author: Elisabeth Beusch
% Compatible with MATLAB R2014a and later
%
% This function computes Halton sequences using the specified prime bases
% and transforms them into standard normal draws. The implementation
% follows chapters 9.3.3–9.3.5 in Train (2003) and includes options for
% scrambling and randomization à la Bhat (2003).
%
% Syntax:
%   [H, Z] = exercisefun(N,dimensions,draws)
%   [H, Z] = exercisefun(N,dimensions,draws,'Name',Value, ...)
%
% Inputs:
%   N           - Number of observational units.
%   dimensions  - Number of dimensions of the integral.
%   draws       - Number of draws per observational unit.
%
% Name-Value Pair Arguments:
%   'prnum'     - Vector of prime numbers used as bases.
%                 Default: ascending primes.
%   'burn'      - Number of initial Halton points to skip.
%                 Default: 50.
%   'leap'      - Number of points to skip between draws.
%                 Default: 0.
%   'random'    - Logical flag to apply randomization à la Bhat (2003).
%                 Set to 1 to enable. Default: 0.
%   'scramble'  - Logical flag to scramble the Halton sequence.
%                 Recommended for high-dimensional settings.
%                 Set to 1 to enable. Default: 0.
%
% Outputs:
%   H - Halton draws of size (N × dimiensions × draws).
%   Z - Corresponding values from a standard normal distribution.
%
% Notes:
%   - The first Halton point (zero) is always dropped.
%   - The function performs basic checks on prime validity and
%     dimensionality.
%
% References:
%   Bhat, C. R., 2003. Simulation estimation of mixed discrete choice
%   models using randomized and scrambled Halton sequences. Transportation
%   Research Part B: Methodological, 37 (9), 837–855.
%
%   Train, K., 2003. Discrete Choice Methods with Simulation. Cambridge
%   University Press.
%
% ---------- BEGIN FUNCTION BODY BELOW ----------

%% Optional input arguments
opts = inputParser;
addParameter(opts,'prnum',0,@isnumeric);
addParameter(opts,'burn',50,@isnumeric);
addParameter(opts,'leap',0,@isnumeric);
addParameter(opts,'random',0,@isnumeric);
addParameter(opts,'scramble',0,@isnumeric);
parse(opts, varargin{:});

prnum      = opts.Results.prnum;
burn       = opts.Results.burn;
leap       = opts.Results.leap;
randhalt   = opts.Results.random;
toscramble = opts.Results.scramble;

%% Define dimensions to be used
if prnum == 0
    p = dimensions;
else
    if size(prnum,1) > size(prnum,2)
        prnum = prnum';
    end

    % Check for consistency and prime validity
    if dimensions ~= size(prnum,2)
        error('Dimensions do not match number of primes supplied');
    end
    if any(~isprime(prnum))
        error('Non-prime was supplied');
    end

    % Determine Halton set dimensionality based on max prime
    p = numel(primes(max(prnum)+1));
end

if dimensions == 1 && (prnum == 0 || prnum == 2)
    error('For dimension == 1, a prime > 2 must be supplied');
end

%% Warnings
if dimensions >= 7 && toscramble == 0
    warning(['With high-dimensional Halton draws, scrambling is ...' ...
        'recommended']);
end

if (burn - max(prnum) <= 10) || (isscalar(prnum) && dimensions >= 13)
    warning('The default burn setting might be too short for your primes');
end

%% Generate Halton sequences
if leap == 0
    h = haltonset(p,'Skip',burn+1); % Sequence starts at zero; skip burn
else
    h = haltonset(p,'Skip',burn+1,'Leap',leap);
end

% Scramble if requested
if toscramble == 1
    h = scramble(h,'RR2');
end

hp = net(h,N*draws);

% Select dimensions based on supplied primes
if prnum ~= 0
    [~, primi] = ismember(prnum,primes(max(prnum)+1));
    hp = hp(:, primi);
end

%% Randomization à la Bhat (Train, 2003, p. 264)
if randhalt == 1
    mu = rand(1,dimensions);
    mu = repmat(mu,N*draws,1);
    hp = hp + mu;
    hp = hp - (hp > 1); % Wrap around values > 1
end

%% Reshape output
H = reshape(hp,draws,N,dimensions); % R x N x J
H = permute(H,[2,3,1]);             % N x J x R

%% Transform to standard normal
if nargout >= 2
    Z = norminv(H);
end