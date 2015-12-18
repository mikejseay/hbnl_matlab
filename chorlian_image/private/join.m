function Y = join(A, B, col)
	if nargin < 3
		col = [1; 1];
	else
		if numel(col) == 1
			col = [col; col];
		else
			col = col(1:2);
		end
	end
	ncA = size(A, 2);
	ncB = size(B, 2);
	Y = [];
	if col(1) > ncA || col(2) > ncB
		return
	end
	A_key = A(:, col(1));
	if length(unique(A_key)) ~= length(A_key)
		return;
	end
	B_key = B(:, col(2));
	if length(unique(B_key)) ~= length(B_key)
		return;
	end	
	A_rest = A(:, setdiff(1:ncA, col(1)));
	B_rest = B(:, setdiff(1:ncB, col(1)));
	[C, AI, BI] = intersect(A_key, B_key);
	Y = [C, A(AI, A_rest), B(BI, B_rest)];



% 	AI = [];
% 	BI = [];
% 	if length(unique(A(:))) ~= length(A(:)) 
% 		fprintf(1, 'argument 1 has duplicate elements\n');
% 		return
% 	end
% 	if length(unique(B(:))) ~= length(B(:))
% 		fprintf(1, 'argument 2 has duplicate elements\n');
% 		return
% 	end	
% 	[C, AI, BI] = intersect(A(:), B(:));
