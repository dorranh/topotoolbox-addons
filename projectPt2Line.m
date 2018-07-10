function [xproj,yproj,dperp,dpar] = projectPt2Line(x,y,linestruct)
    % Be sure to avoid any trailing nan values that ArcGIS tends to create
    if isnan(linestruct.Y(end)) || isnan(linestruct.X(end))
        pt1 = [linestruct.X(1);linestruct.Y(1)];
        pt2 = [linestruct.X(end-1);linestruct.Y(end-1)];
    else
        pt1 = [linestruct.X(1);linestruct.Y(1)];
        pt2 = [linestruct.X(end);linestruct.Y(end)]; 
    end
    
    % 1.) Get 2 Points in this line
    slp = (pt2(2) - pt1(2)) /(pt2(1) - pt1(1));
    shft = [1;slp];
    % 2.) Project basis vectors [1,0],[0,1] to translated line
    P_2 = zeros(2,2);
    P_2(:,1) =  (dot([1;0],shft)/dot(shft,shft))*shft;
    P_2(:,2) =  (dot([0;1],shft)/dot(shft,shft))*shft;
    % 3.) Construct Projection Matrix, P
    P = zeros(3,3);
    P(1:2,1:2) = P_2;
    P(3,:) = [0 0 1];
    P(1:2,3) = pt1 - P_2*pt1;
    
    % Get the xy coord of this pt in the reference  line
    projected_xy = P*[x;y;1];
    xproj = projected_xy(1);
    yproj = projected_xy(2);
    % Get the vector from this point to the reference line
    newxydir = [projected_xy(1)-x,projected_xy(2)-y];
    % Determine it's magnitude (their distance)
    dperp = norm(newxydir);
    % Get distance from beginning of line to projected point (parallel
    % dist)
    dpar = sqrt((yproj-pt1(2))^2 + (xproj-pt1(1))^2); 
end