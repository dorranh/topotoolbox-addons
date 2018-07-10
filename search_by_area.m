function [outlets, areas] = search_by_area(search_init, min_area, max_area, DEM, varargin)
%search_by_area finds drainage basins with areas that fall within the range
%of inputs 'min_area' and 'max_area' given a starting index in your DEM,
%'search_init'.
%
%   This function performs a recursive search for drainage basins for each initialization
%   point by walking up the catchment, following tributaries until either a
%   drainage basin with an acceptable area is found or the channel head is
%   reached.  Output outlets will be located one pixel upstream from a
%   confluence.
%
%-- Arguments: --
%
% * Required *
%   search_init    : a single value or array of linear indices into the DEM
%                    where you want to start a search. 
%   min_area       : minimum acceptable drainage area for a catchment
%   max_area       : maximum acceptable drainage area for a catchment
%   DEM            : GRIDobj instance of a DEM with sinks filled
%   
% * Optional *     [Note: 'stream' and 'area' must be provided together]
%   stream         : STREAMobj instance to use for search
%   area           : GRIDobj instance of a drainage area grid
%   minarea_stream : minimum drainage area for defining a STREAMobj 
%
%-- Returns: --
%
%   outlets : cell array of length(search_init), with each element
%             containing an array of corresponding outlet linear indices
%   areas   : cell array of length(search_init), with each element
%             containing an array of corresponding drainage basin areas
%
%-- Usage: --
%
%   1.) Pass in a DEM and the function will compute a Flow
%       Accumulation Grid, FLOWobj, and STREAMobj Automatically:
%       >> [outlets,areas] = search_by_area(search_init, min_area, max_area, DEM);
%
%   OR
%
%   2.) Pass in a pre-computed drainage area grid (can be created using the output
%       from flowacc(...) and stream object:
%       >> [outlets,areas] = search_by_area(search_init, min_area, max_area, ...
%                                           DEM, 'stream', S, 'area', A);
%
%-- Dependencies: --
% TopoToolbox (Developed for Version 2.2)
%
% Written By: Dorran Howell (dorran.howell@gmail.com)
% Last Updated 03/14/2018

p = inputParser;
p.FunctionName = 'search_by_area';
% Required Inputs
addRequired(p,'search_init', @(x) isscalar(x) | isvector(x));
addRequired(p,'min_area', @(x) isscalar(x));
addRequired(p,'max_area', @(x) isscalar(x));
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
% Optional Inputs
addOptional(p,'stream', false, @(x) isa(x,'STREAMobj'));
addOptional(p,'area', false, @(x) isa(x,'GRIDobj'));
addOptional(p,'minarea_stream', 1e6, @(x) isscalar(x));
%addOptional(p,'parallel',false, @(x) islogical(x));

% Parse Inputs
parse(p,search_init,min_area,max_area,DEM,varargin{:});
S = p.Results.stream;
A = p.Results.area;
minarea_stream = p.Results.minarea_stream;

% Prepare area grid and stream object
if ~isa(A,'GRIDobj') || ~isa(S,'STREAMobj')
    fprintf('No Area Raster or STREAMobj Provided. If your DEM is large it is\nrecommended that you precalculate these.\nComputing a STREAMobj, Flow Accumulation, and Drainage Area Grids...\n')
    FD = FLOWobj(DEM);
    FA = flowacc(FD);
    A = FA .* DEM.cellsize^2.0;                    
    S = STREAMobj(FD,'minarea',minarea_stream/DEM.cellsize^2);
    fprintf('Done.\n')
end

% Get Locations of Confluences
fprintf('Identifying Confluences...\n');
confluences = streampoi(S,'confluences','ix');
fprintf('Done.\n');

% Prepare STREAMobj vectors for search function
% Convert S.ix and S.ixc to be in terms of linear indices into DEM
S_dem_ix = S.IXgrid(S.ix);
S_dem_ixc = S.IXgrid(S.ixc);
S_dem_ix = flip(S_dem_ix);
S_dem_ixc = flip(S_dem_ixc);
IXgrid_bottomup = flip(S.IXgrid);

% Run recursive search for each initialization point and append results to
% output
fprintf('Starting Drainage Basin Search...\n')
outlets = cell([length(search_init),1]);
areas = cell([length(search_init),1]);
total_found = 0;
% Snap points to stream
[init_x,init_y] = ind2coord(DEM,search_init);
[~,~,snapped] = snap2stream(S,init_x,init_y);

update_indices = 1:round(length(search_init)/100):length(search_init);
for i = 1:length(search_init)
    if ismember(i,update_indices)
       fprintf("Searching Outlet %d of %d...\n",i,length(search_init)); 
    end
    this_init_point = snapped(i);
    [ix,db_area] = check_upstream(this_init_point,S_dem_ix,S_dem_ixc,IXgrid_bottomup,A,confluences,min_area,max_area);
    outlets{i} = ix;
    areas{i} = db_area;
    total_found = total_found + length(ix);
end
fprintf('Search Complete. Found a total of %d drainage basins for all initialization points.\n',total_found)

end
