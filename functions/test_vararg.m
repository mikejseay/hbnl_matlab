function varargout = test_vararg(varargin)

% this function tries to read the name of the varargout... but it can't be
% done

in1=varargin{1};

if isscalar(in1)
    varargout{1}=2*in1;
end

end
