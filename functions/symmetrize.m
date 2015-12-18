function v_out = symmetrize(v_in, center)
% make a 2-element vector (of plotting limits) symmetric about its center

if nargin<2
    if sign(v_in(1)) ~= sign(v_in(2))
        center=0;
    else
        center=mean(v_in);
    end
end

if abs(v_in(1) - center) > abs(v_in(2) - center)
    v_out(1) = v_in(1);
    v_out(2) = center + abs(v_in(1) - center);
else
    v_out(2) = v_in(2);
    v_out(1) = center - abs(v_in(2) - center);
end

end