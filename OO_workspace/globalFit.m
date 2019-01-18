function [pfit,err,gfOpts] = globalFit(dataMatrix,w1,w3,p_array,gfOpts,lb,ub)
% [x,xError,struct] = globalFit(dataMat,w1,w3,p0,gfOptions

% Scale intensity of the calculated Response Function to match that of the
% 2D-IR spectrum. This also adds some natural weighting (by the signal size
% of each spectrum) of the fit.
scale = gfScale(dataMatrix);

% The variable MAXSCALE allows us to renormalize the data, so that we aren't constantly
% exceeding our tolerance because of the large amplitudes of our spectra.
maxScale = abs(min(min(min(dataMatrix))));

err_fun = @(p,gfOpts) sum(sum(sum((dataMatrix-analyticalResponseFunctionsFun(p,w1,w3,gfOpts).*scale).^2)))./maxScale;

%tolfun relates to the squared residual in the error
%tolx relates to the displacement in your parameters
tolfun = 1e-17;
tolx = 1e-10;
maxfun = 1e3; %maximum number of function evaluations


% set fit parameters (versions of MatLab before 2013a did not have the OPTIMOPTIONS function)
if exist('optimoptions')==2,
    opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
elseif exist('optimset')==2,
%    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
else
    warning('Could not set parameters! Look for the right options functon for your installation.');
end

opt = optimoptions(opt,'UseParallel',true);

tic

% fmincon has required parameters of error function and initial guess. We
% are omitting any linear inequalities that we could subject the minimizer
% to (e.g. we could include A and b (input parameters 3 & 4), and would
% require that A*pfit <= b.
[pfit,err] = fmincon(err_fun,p_array,[],[],[],[],lb,ub,[],opt,gfOpts)

toc

gfOpts.pfit = pfit;
gfOpts.fitMethod = 'unweighted fit';

for ii = 1:length(pfit)
    fprintf(1,'%20s\t%12f\n',gfOpts.pnames{ii},pfit(ii));
end