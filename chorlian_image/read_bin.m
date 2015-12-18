function [A, count] = read_bin(varargin)
	file =  varargin{1};
	fid = fopen(file, 'r', 'b');
	if fid == -1
		A = [];
		count = 0;
	else
		if nargin == 3
			nr = varargin{2}; nc = varargin{3};
			[A, count] = fread(fid, [nr, nc], 'float32');
		else
			[A, count] = fread(fid, inf, 'float32');
		end
		fclose(fid);
	end
	
	
% function A = read_bin(file)
% 	fid = fopen(file, 'r');
% 	if fid == -1
% 		A = [];
% 	else
% 		A = fread(fid, inf, 'float32');
% 		fclose(fid);
% 	end
		
	
