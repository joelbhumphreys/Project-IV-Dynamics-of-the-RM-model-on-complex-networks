function [P, t, num_crossings] = poincare_map(X, plane)
    % P = poincare_map( X [, plane] )
    % given N-dimensional time-series data X,  
    % find points of the time series that cross a given plane. 
    %   X ( t, variables ) is a T x N matrix of N-dimensional states evolving over time
    %   plane.norm = N-dim normal vector  (default [1,0,0,0...]
    %   plane.dist = distance from origin (default 0)
    % sanjay g manohar 2019

    FAST_MODE = false; % calculate in matrix form, but do not interpolate points?

    % define the poincare plane
    if ~exist('plane','var')
        plane.norm = zeros(size(X,2),1);  % column vec
        plane.norm(1) = 1;
        plane.dist = 0;
    end
    plane.norm = plane.norm/norm(plane.norm); % ensure length = 1

    % distance of each state point from the plane
    % project the state onto the plane's normal vector (dot product)
    projections = X * plane.norm - plane.dist;
    if FAST_MODE 
        % select points where the plane is crossed
        crossings = projections(1:end-1) < 0 & projections(2:end) > 0 ;
        P = X(crossings,:); 
        t = find(crossings); 
        num_crossings = sum(crossings);
    else % SLOW MODE: estimate each plane-crossing using linear interpolation
        P = [] ; % keep a list of plane-crossing points
        t = [] ; % keep a list of the times of crossing
        num_crossings = 0;
        prev_crossings = []; % store the previously added crossing points
        for i=2:size(X,1) % for each timepoint
            p1 = projections(i-1); % previous state's distance from plane
            p2 = projections(i);   % current state's distance from plane
            % is this a plane-crossing?
            if p2>0 && p1<0
                % linear interpolation to mix the two points, to estimate the crossing
                % point
                mix = -p1/(p2-p1); % how much of p2 to mix into p1?
                Xp = X(i-1,:) + mix * (X(i,:)-X(i-1,:)); % estimated plane crossing
                
                % check if the new crossing point is sufficiently far away
                % from every previous crossing point
                if isempty(prev_crossings) || all(min(pdist2(Xp, prev_crossings)) > 1e-2)
                    % add a row to P
                    P = [P; Xp];
                    t = [t; i];
                    num_crossings = num_crossings + 1;
                    prev_crossings = [prev_crossings; Xp]; % update the previous crossing points
                end
            end
        end
    end
end
