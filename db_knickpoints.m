function [basinShape, streamShape, DEMc] = db_knickpoints(DEM,FD,ix,dbid,varargin)
% basinStatisticsFunction calculates key statistics for drianage basins
% for user defined drainage basins the the input DEM
% Inputs:
%   DEM: [req] Grid containing elevation data for an area [TopoToolbox.GRIDObj]
%   ix: [req] int or int array of DB outlet linear indices for the input DEM
%   mn: [req] m/n ratio for Chi calculation (see Seans_functions.ChiFits)
%   critA: [req] Critical drainage area for Chi calculation 
%   Ao: [req] Reference area for Chi calculation
%   hypsBin: [opt,default=50] Binning increment [m] for Basin Hypsometry
%
% Outputs:
% dataTable - Table containing output statistics for each drainage basin
%
% Written By: Dorran Howell
% ________________________________________________________________________

%% Parse Arguments
p = inputParser;
p.FunctionName = 'db_knickpoints';
% Args
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
addRequired(p,'ix', @(x) isscalar(x) || isvector(x));
addRequired(p,'dbid',  @(x) isscalar(x) || isvector(x));
% Kwargs
addParameter(p,'tol',25,@(x) isscalar(x));
addParameter(p,'critA',1e6,@(x) isscalar(x));
addParameter(p,'critLength',3e3,@(x) isscalar(x))
addParameter(p,'plot',false,@(x) islogical(x))

parse(p,DEM,FD,ix,dbid,varargin{:});

%% Pre-processing
% Initialize temporary array for containing statistics
data = struct();

for i=1:length(ix)
% Get mask for this DB
DEMc = DEM;
DB = drainagebasins(FD,ix(i));
maskc = crop(DB,DB);
DEMc = crop(DEMc,DB);
% Create drainage clipped flow objects
FDdb = FLOWobj(DEMc);
Adb = flowacc(FDdb);
% Transform linear indices into the clipped DEM frame of reference
[ox,oy] = ind2coord(DEM,ix(i));
ixdb = coord2ind(DEMc,ox,oy);
% Generate a new stream object for the clipped area
Sdb = STREAMobj(FDdb,'unit','mapunits','minarea',p.Results.critA,'outlets',ixdb);
% Null out areas in DEM outside of drainage basin mask
DEMc.Z(~maskc.Z) = nan;
% Ensure that the cropped flow accumulation and elevation grids are
% alligned
validatealignment(Adb,DEMc);
% Extract main drainage network for this DEM 
Sdb = klargestconncomps(Sdb);
% Remove short stream segments
Sdb = removeshortstreams(Sdb,p.Results.critLength);
% Run carving to reduce noise in profile
zs = quantcarve(Sdb,DEMc,.5,'split',false);
% Search for knickpoints
[zk,kp] = knickpointfinder(Sdb,zs,'split',false,'plot',false,'tol',p.Results.tol,'verbose',false);

% order =  Sdb.orderednanlist;
% d     = nan(size(order));
% I     = ~isnan(order);
% d(I)  = dist(order(I));
if p.Results.plot
    figure(1)
    clf
    scatter(distance(Sdb),zs)
    hold on
    scatter(kp.distance,kp.z,kp.dz,'sk','MarkerFaceColor','r')
    hold off
end

end
