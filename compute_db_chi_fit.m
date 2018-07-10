% E.g. [chi,r2,c,z] = compute_db_chi(DEM,FD,outlets,0.45,1e6)

function [chi_grid, r2,z,chi_vec] = compute_db_chi_fit( DEM, A, S, outlet_ix, mn, Ao,trim)

% Make non DB areas nan if trim is TRUE
if trim
    % Get the outline indices of this drainage basin
    DB = drainagebasins(FD,outlet_ix);

    % Copy elevation data and trim to drainage basin extent
    DEM.Z(DB.Z < 1) = nan;
end

% Create CHI Map for this drainage basin
% I am a lazy bastard.  Stolen from Sean's chi_profiler functions

% Calculate Drainage Area Raster
cs = DEM.cellsize;
AREA   = A.*(cs^2);

% ----------BEGIN SEANS CHI MAP STUFF ------------
% Declare STREAMobj variables for faster processing through forloop
ordList = S.orderednanlist;
strmBreaks = find(isnan(ordList));
GridID = S.IXgrid;

Sz = double(DEM.Z(GridID));         % elevation

% get variables ready for chi integration
chis = zeros(size(S.distance));
Six = S.ix;
Sixc = S.ixc;
Sx = S.distance;
Sa = (Ao./(AREA.Z(S.IXgrid))).^mn;

% calculating chi for the entire river network
for lp = numel(Six):-1:1
    chis(Six(lp)) = chis(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
end

ChiGrid = DEM;
ChiGrid.Z = nan(size(DEM.Z));
ChiGrid.Z(S.IXgrid) = chis;

chi_grid = ChiGrid;
S_db = S;

chi_vec = ChiGrid.Z(S.IXgrid);
z = DEM.Z(S.IXgrid);
res = regstats(z,chi_vec,'linear','rsquare');
r2 = res.rsquare;

end

