function Chi = PWeightedChiFits(S,DEM,A,P,varargin)
%
% ChiFits.m will find the best-fit m/n (concavity) and corresponding ks
% (channel steepness index) using the intergal approach (Perron and Royden,
% 2013). The best-fit m/n is found two different ways: (1) maximizing the 
% r-squared of a linear fit to basin wide chi-elevation data, and (2)
% minimizing the variance in chi-elevation data based on 100 evenly space
% bins in chi-space.
%
% ChiFits.m also calculates chi for the river network based on the two
% best-bit m/n values. It will calculate the drainage basin wide 
% normalized steepness index (ksn) and chi provided a reference m/n value.
%
% ChiFits.m is basically a modified version of the TopoToolbox function
% chiplot.m by Wolfgang Schwanghart (Schwanghart and Scherler, 2013).
%
% Inputs:
% Required:
% 1) S --> STREAMobj. The stream network must be consist of only one
%          connected component (e.g. it must only have one outlet
% 2) DEM --> GRIDobj
% 3) A --> flow accumulation calculated by flowacc (GRIDobj)
%
% Optional parameter name/value pairs {default}:
% 1) 'a0' --> reference drainage area for chi integration {1}
% 2) 'mn' --> reference m/n (concavity) {0.5}
% 3) 'mnplot' --> plot chi-elevation for mn [0.1:0.1:0.9]? {true}, false
% 4) 'plot' --> plot mn residuals for [0.001:0.001:0.999]? {false}, true
%
% Outputs:
% 1) Chi --> structure array that contains...
%       .mn_r2: best-fit m/n based on maximum r-squared.
%       .mn_var: best-fit m.n based on minimum variance.
%       .ks_r2: steepness index from .mn_r2.
%       .ks_95.uc_r2: 95% confidence interval for .ks_r2.
%       .ks_var: steepness index from .mn_var.
%       .ks_95.uc_var: 95% confidence interval for .ks_var.
%       .mnref: reference m/n (convavity)
%       .ksn: normalized steepness index from .mnref.
%       .ksn_95uc: 95% confidence interval for .ksn.
%       .chi_r2: vector of chi calculated from .mn_r2.
%       .chi_var: vector of chi calculated from .mn_var.
%       .chi_mnref: vector of chi calculated from .mn_mnref.
%       .z: vector of channel elevation relative to mean sea level.
%       .zb: vector of channel elevation relative to outlet.
%       .x: x-coordinates of network
%       .y: y-coordinates of network
%       .distance: vector of channel distances
%       .IXgrid: linear indicies of river network in DEM GRIDobj
%       .A: network vector of drainage area in units of cells
%
% Example
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S);
%     C   = ChiFits(S,DEM,flowacc(FD));
%
% Author: Sean F. Gallen
% Date modified: 02/14/2017 <3<3<3<3<3<3<3<3
% email: sean.gallen[at]erdw.ethz.ch


% Parse Inputs
p = inputParser;         
p.FunctionName = 'ChiFits';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'A', @(x) isa(x,'GRIDobj'))
addRequired(p,'P', @(x) isa(x,'GRIDobj'))


%addParameter(p,'mn',[],@(x) isscalar(x) || isempty(x));
addOptional(p,'plot',true,@(x) isscalar(x));
addOptional(p,'mnplot',false,@(x) isscalar(x));
addOptional(p,'a0',1,@(x) isscalar(x) && isnumeric(x));
addOptional(p,'p0',1,@(x) isscalar(x) && isnumeric(x));
addOptional(p,'mn',[],@(x) isscalar(x) || isempty(x));

parse(p,S,DEM,A,P,varargin{:});
S     = p.Results.S;
DEM   = p.Results.DEM;
A     = p.Results.A;
a0    = p.Results.a0;
p0    = p.Results.p0;
mnplot = p.Results.mnplot;
to_plot = p.Results.plot;

outlets = streampoi(S,'outlet','ix');

if length(outlets) > 1
    disp('Note: Your STREAMobj has more than 1 outlet...');
    disp('Only using largest connected stream network in basin');
    S = klargestconncomps(S,1);
    outlets = streampoi(S,'outlet','ix');
end


% get elevations
zx = double(DEM.Z(S.IXgrid));
zb = double(DEM.Z(outlets));

% NEW - get precipitation weighed drainage area
a = double(   (a0*p0) ./ ((A.Z(S.IXgrid)*(A.cellsize^2)) .* P.Z(S.IXgrid))   );

% get distance
x = S.distance;

% loop through and calculate chi based on different mn values.
% calculate mean variance of the network and the r-squared based on a
% linear regression through chi-z data.
mntest = 0.001:0.001:0.999;

% allocate memory for variables
chitest = zeros(length(x),length(mntest));
mnVar = nan(length(mntest),1);
mnR2 = nan(length(mntest),1);

for k = 1:length(mntest);
    % calculate chi
    chitest(:,k) = cumtrapz(S,a.^mntest(k));
    chiTemp = chitest(:,k);
    
    % calculate number of bins
    bins = linspace(min(chiTemp),max(chiTemp),100);
    
    % allocate memory for variance vector
    chivar= nan(length(bins),1);
    for j = 1:length(bins)-1;
        chivar(j) = nanvar(zx(chiTemp >= bins(j) & chiTemp < bins(j+1)));
    end
    mnVar(k) = nanmean(chivar(chivar > 0));
    
    CC = [ones(size(zx)) chiTemp];
    [b,bint,r,rint,stats] = regress(zx-zb,CC);
    mnR2(k) = stats(1);
end

if to_plot
    figure()
    subplot(1,2,1)
    plot(mntest,mnR2,'b-','lineWidth',1.5); hold on
    plot(mntest(mnR2 == max(mnR2)),mnR2(mnR2 == max(mnR2)),'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8]);
    xlabel('m/n'), ylabel('r-squared');
    title('R-squared from linear regression')
    %legend({'results',['best-fit m/n: ',num2str(mntest(mnR2 == max(mnR2)))]});


    subplot(1,2,2)
    plot(mntest,mnVar,'b-','lineWidth',1.5); hold on
    plot(mntest(mnVar == min(mnVar)),mnVar(mnVar  == min(mnVar)),'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8]);
    xlabel('m/n'), ylabel('mean variance');
    title('Mean variance of elevation per \chi bin');
    %legend({'results',['best-fit m/n: ',num2str(mntest(mnVar == min(mnVar)))]});
end

% plot results every 0.1 mn
if mnplot
    figure()
    plot(chitest(:,rem(mntest,0.1) == 0),zx-zb,'.');
    xlabel('\chi [m]')
    ylabel('elevation [m]');
    title('\chi plots for different values of mn')
end

Chi.mn_r2 = mntest(mnR2 == max(mnR2));
Chi.mn_var = mntest(mnVar == min(mnVar));
Chi.mn_r2 = Chi.mn_r2(1);
Chi.mn_var = Chi.mn_var(1);
chi_r2 = cumtrapz(S,a.^Chi.mn_r2);
chi_var = cumtrapz(S,a.^Chi.mn_var);
    
CC = [ones(size(zx)) chi_r2];
[b,bint,k,rint,stats] = regress(zx-zb,CC,0.05);
Chi.ks_r2 = b(2)*a0^Chi.mn_r2;
Chi.ks_95uc_r2 = (bint(4)-bint(2))/2;

CC = [ones(size(zx)) chi_var];
[b,bint,k,rint,stats] = regress(zx-zb,CC,0.05);
Chi.ks_var = b(2)*a0^Chi.mn_var;
Chi.ks_95uc_var = (bint(4)-bint(2))/2;

% calculate a ksn value
if isempty(p.Results.mn);
    mnref = 0.5;
else
    mnref = p.Results.mn;
end

chi_mnref = cumtrapz(S,a.^mnref);
CC = [ones(size(zx)) chi_mnref];
[b,bint,k,rint,stats] = regress(zx-zb,CC,0.05);
Chi.mnref = mnref;
Chi.ksn = b(2)*a0^mnref;
Chi.ksn_95uc = (bint(4)-bint(2))/2;

Chi.chi_r2 = chi_r2;
Chi.chi_var = chi_var;
Chi.chi_mnref = chi_mnref;
Chi.z = zx;
Chi.zb = zx-zb;
Chi.x = S.x;
Chi.y = S.y;
Chi.distance = S.distance;
Chi.IXgrid = S.IXgrid;
Chi.A = A.Z(S.IXgrid);
end

    
