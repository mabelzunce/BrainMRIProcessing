function [h,p,test,df] = chi2test2(x1,x2,alpha)
%CHI2TEST2  Two-sample chi-squared test for discrete data.
%   H = CHI2TEST2(X) performs a chi-squared two-sample test of the
%   similarity of two samples, under the null hypothesis that two samples
%   come from the same (non-specified) distribution.
%
%   X1 and X2 are two vectors containing discrete (categorical) data.
%   CHI2TEST2 treats all unique values as separate 'bins' for testing
%   purposes. X1 and X2 can be different lengths. CHI2TEST2 treats NaNs as
%   missing values, and ignores them.
%
%   H = TTEST(...,ALPHA) performs the test at the significance level
%   (100*ALPHA)%.  ALPHA must be a scalar. Default value of ALPHA is 0.05.
%   i.e., H = 1 suggests different distributions.
%
%   [H,P] = TTEST(...) returns the p-value, i.e., the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis is true.  Small values of P cast doubt on the validity of
%   the null hypothesis (i.e. suggests different distributions).
%
%   [H,P,TEST] = TTEST(...) returns TEST, the value of the test statistic.
%
%   [H,P,TEST,DF] = TTEST(...) returns DF, the degrees of freedom of the test.
%
%   EXAMPLE 1: Generate two uniform random vectors, should be H = 0
%       rng('default')
%       x1 = round(rand(1000,1)*10);
%       x2 = round(rand(1000,1)*10);
%       [h,p,test,df] = chi2test2(x1,x2)
%             h =
%                  0
%             p =
%                 0.4505
%             test =
%                 9.8868
%             df =
%                 10
% 
%   EXAMPLE 2: Generate two random vectors, one uniform and one normal,
%   should be H = 1
%       rng('default')
%       x1 = round(rand(1000,1)*10);
%       x2 = 5 + round(randn(1000,1));
%       [h,p,test,df] = chi2test2(x1,x2)
%             h =
%                 1
%             p =
%                 0
%             test =
%                 758.0990
%             df =
%                 10
%                 
%
%   Reference:
%   http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/chi2samp.htm
%           Accessed March 29, 2013, and cites: "Numerical Recipes in
%           Fortan: The Art of Scientific Computing", Second Edition,
%           Press, Teukolsky, Vetterlling, and Flannery, Cambridge
%           University Press, 1992, pp. 614-622. 
%
%   Copyright 2013 James Meldrum. Provided for entertainment only; author
%   is not responsible for any use or misuse of code. Verify before use!
%   Revision: 1, 2013/03/29 $

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha) || alpha <= 0 || alpha >= 1
    error('Bad alpha value, please try again.');
end

if isempty(x1) || isempty(x2)
    error('Missing or empty X1 or X2 vectors, please try again.');
end

% Purge input vectors of any NaN
x1 = x1(~isnan(x1)); 
x2 = x2(~isnan(x2)); 

% Calculate frequencies and scaling constants
bins = unique([x1(:,1); x2(:,1)]); % create a bin for each unique value
freq1 = histc(x1,bins);
freq2 = histc(x2,bins);
k1 = sqrt(sum(freq2)/sum(freq1));
k2 = sqrt(sum(freq1)/sum(freq2));


% Calculate test statistic
test = sum((k1*freq1-k2*freq2).^2./(freq1+freq2));

% Calculate degrees of freedom
if size(x1,1)==size(x2,1)
    df = length(bins)-1;
else
    df = length(bins);
end

% Calculate p-value under chi-squared distribution
p = 1 - chi2cdf(test,df);

% Test calculated versus desired significance level
h = cast(p <= alpha, class(p));
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha
