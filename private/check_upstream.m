function [ix_db,area_db] = check_upstream(n,ix,ixc,ixgrid,A,confluences,minArea,maxArea)
    ix_db = [];
    area_db = [];
    start = find(ixc==n);
    this_point = n;
    % Work your way upstream
    for i = start:length(ixc) 
        % If this is a confluence, get the adjacent upstream points in all
        % tributaries
        if any(confluences == this_point)
           forks = find(ixc == this_point); 
           forks_ix = ix(forks);
    
           % Gather areas for the node upstream from the confluence on each
           % tributary
           A_forks = A.Z(forks_ix);
           test = A_forks >= minArea & A_forks <= maxArea;
           too_big = A_forks > maxArea;
           
           % Append any outlets that fall within our range
           if any(test)
               ix_db = [ix_db;forks_ix(test)];
               area_db = [area_db;A_forks(test)];
           end
           
           % If there are no more DBs that can be looked at quit
           if ~any(too_big)
                break
           else
               % Otherwise keep search in the remainings DBs
               remaining_dbs = forks_ix(too_big);
               for j = 1:length(remaining_dbs)
                   this_n =remaining_dbs(j);
                   % Make recursive call
                   [sub_db_ix,sub_db_A] = check_upstream(this_n,ix,ixc,ixgrid,A,confluences,minArea,maxArea);
                   % Append outputs
                   ix_db = [ix_db;sub_db_ix];
                   area_db = [area_db;sub_db_A];
               end
               break
           end
           break
        %If this is not a confluence, move to the next point upstream and
        %repeat
        else
            next_point = ix(ixc==this_point);
            % End once you hit the channel head
            if isempty(next_point)
               break 
            end
            this_point = next_point;
        end
    end   
end
    