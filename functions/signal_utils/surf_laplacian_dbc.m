function lap = surf_laplacian_dbc(data, coords, radius, n_points)
	%  data is [n_times X n_chans]
	%  coords is [n_chans X 3]  which are unit vectors
	%	radius and n_points are scalars
	%	if radius is non-zero, use all points within radius
	%	with minimum of n_points otherwise
	%  n_points determines number of nearest neighbors to use  
	%  output is the laplacian at the coords
	%	note we could calculate Laplacian at different points
	%  also, that we could compute interpolated values of potential
	%  this would require slight modifications of the 
	%	find_nearest_neighbors() and get_local_coords() functions
	%  method is to compute nearest neighbors t0 electrode
	%  make local coordinate system for that electrode
	%	calculate design matrix for that local coordinate system
	%  solve for all time points simultaneously
	%  Copyright (C) 2005 by David B. Chorlian
	%%% Written by David B. Chorlian, based on the method described in: Wang K, Begleiter H (1999) Local polynomial estimate of surface Laplacian. Brain Topogr 12:19-29. PMID: 10582562.


% 	This program is free software; you can redistribute it and/or
% 	modify it under the terms of the GNU General Public License as
% 	published by the Free Software Foundation; either version 2 of
% 	the License, or (at your option) any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% 	GNU General Public License for more details.
% 
% 	You should have received a copy of the GNU General Public
% 	License along with this program; if not, write to the Free
% 	Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
% 	Boston, MA  02110-1301, USA

% 		
% 			Problem: Given data y, and geometrically determined (by
% 			location of point to be fit) weights w and design matrix A, 
% 			find coeffients x such that
% 			w * A * x = w * y
% 			
% 			A is determined by the number of data points to be used in
% 			the fit, and by their location, given that a local quadratic
%			is to be used for the approximation
% 			w is determined by the distance of the data points from the
% 			point to be fit.

% 			use Matlab backslash divide to solve w * A * x = w * y
%			From Matlab documentation:
% 			If A is not square, then Householder reflections
% 			are used to compute an orthogonal-triangular
% 			factorization. A*P = Q*R where P is a permutation,
% 			Q is orthogonal and R is upper triangular (see
% 			qr). The least squares solution X is computed with
% 			X = P*(R\(Q'*B)) MATLAB uses the LAPACK routines.
%			(This last is easy because R is upper triangular)

% 			
% 			1) obtain A and w from local coordinates
% 			2) calculate w * A
% 			4) calculate w * y
% 			5) calculate x = mldivide(w * A, w * y)
	
	[n_times, n_chans] = size(data);
	% coord_dotprod used in subsequent calculation of local coordinates
	coord_dotprod = coords * coords';  
	coord_distances = 1 - coord_dotprod;
	lap = zeros(n_times, n_chans);
	for m = 1:n_chans
		nearest_neighbors = ...
			find_nearest_neighbors(coord_distances(:, m), ...
			radius, n_points);
		[local_coords, weights] = ...
			get_local_coords(coords, coord_dotprod, coord_distances, ...
				nearest_neighbors);
		wA = make_design_matrix(local_coords, weights);
		lap(:, m) = ...
			laplacian(data, nearest_neighbors(1:(end - 1)), wA, weights);
	end


function [nearest_neighbors, radial_distance] = ...
	find_nearest_neighbors(distance, radius, n_points)
	% distance obtained from 1 - coords * coords';
	[x, Ix] = sort(distance);
	if radius > 0
		x_len = length(find(x < radius));
		if x_len > n_points
			n_points = x_len - 1;
		end
	end
	% 1st point is center 
	% center + n neighbors + 1 to calculate weights
	radial_distance = x(n_points + 1);
	nearest_neighbors = Ix(1:(n_points + 2));
	
function [local_coords, weights] = ...
	get_local_coords(coords, dotprods, distance, nearest_neighbors)
	% center point is coord(nearest_neighbors(1))
	% construct local x, y system in tangent plane to sphere at center point
	n_points = length(nearest_neighbors) - 1;
	local_coords = zeros(n_points, 2);
	weights = ones(n_points, 1);
	local_coords(1, :) = 0;
	m = nearest_neighbors(1);

%	calculate projection of coord(n) on tangent plane, 
% 		the projection of pt_n onto the x-y plane is 
% 			proj_pt_n = pt_n - <ctr_pt, pt_n> * ctr_pt
% 		so the local coordinates are
% 			proj_pt_n.x  = <proj_pt_n, local_x_axis>
% 			proj_pt_n.y  = <proj_pt_n, local_y_axis>	


%  account for special cases of coordinates of center point
%  following case is unnecessary; accounted for by default	
% 	if (abs(coords(m, 3)) + abs(coords(m, 1))) < eps % coord is 0, 1, 0
% 		local_x_axis = [1, 0, 0];
% 		local_y_axis = [0, 0, 1];
	if (abs(coords(m, 2)) + abs(coords(m, 1))) < eps % coord is 0, 0, 1
		local_x_axis = [1, 0, 0];
		local_y_axis = [0, 1, 0];
	elseif (abs(coords(m, 3)) + abs(coords(m, 2))) < eps % coord is 1, 0, 0
		local_x_axis = [0, 1, 0];
		local_y_axis = [0, 0, 1];	
	elseif abs(coords(m, 1)) < eps % coord is 0, y, z
		local_x_axis = [1, 0, 0];
		local_y_axis = [0, coords(m, 3),  -coords(m, 2)];
%  following case is unnecessary; accounted for by default
% 	elseif abs(coords(m, 2)) < eps % coord is x, 0, z
% 		local_x_axis = [0, 1, 0];
% 		local_y_axis = [coords(m, 3), 0,  -coords(m, 1)];
	elseif abs(coords(m, 3))  < eps % coord is x, y, 0
		local_x_axis = [0, 0, 1];
		local_y_axis = [coords(m, 2), -coords(m, 1), 0];	
	else
		local_x_axis = [coords(m, 3), 0, -coords(m, 1)];
		local_y_axis = [-coords(m, 1) * coords(m, 2), ...
			coords(m, 1)^2 + coords(m, 3)^2, -coords(m, 3) * coords(m, 2)];
	end
	local_y_axis = local_y_axis'/(norm(local_y_axis));
	local_x_axis = local_x_axis'/(norm(local_x_axis));
	% calculate projection of coord(n) on tangent plane, 
% 		the projection of pt_n onto the x-y plane is 
% 			proj_pt_n = pt_n - <ctr_pt, pt_n> * ctr_pt
% 		so the local coordinates are
% 			proj_pt_n.x  = <proj_pt_n, local_x_axis>
% 			proj_pt_n.y  = <proj_pt_n, local_y_axis>	
% note that there are different choices for calculating
% weight vector
	
	extremal_distance = ...
		(distance(m, nearest_neighbors(n_points)) + ...
			distance(m, nearest_neighbors(n_points + 1)))/2;
	for n = 2:n_points
		mm = nearest_neighbors(n);
		proj = coords(mm, :) - dotprods(m, mm) .* coords(m, :);
		local_coords(n, 1) = proj * local_x_axis;
		local_coords(n, 2) = proj * local_y_axis;
		weights(n) = (1 - (distance(m, mm)/extremal_distance)^2)^2;
	end
	
function wA = make_design_matrix(local_coords, weights)
	n_points = size(local_coords, 1);
	A = ones(n_points, 6);
	A(:, 2) = local_coords(:, 1);
	A(:, 3) = local_coords(:, 2);
	A(:, 4) = local_coords(:, 1) .^ 2;
	A(:, 5) = local_coords(:, 2) .^ 2;
	A(:, 6) = local_coords(:, 1) .* local_coords(:, 2);
	wA = A .* repmat(weights, [1, 6]);
	
	
function [lX, X] = laplacian(data, nearest_neighbors, wA, weights)
	% data is [n_samps, n_chans]
	% wA is [n_chans, 6]
	% weights is [n_chans, 1]
	% we want to solve wA * x = w * data
	% where x is the vector of coefficients of the design matrix
	n_samps = size(data, 1);
	U = (data(:, nearest_neighbors) .* repmat(weights', [n_samps, 1]))';
	X = mldivide(wA, U);
	lX = -((X(4, :) + X(5, :)))' * 2;
