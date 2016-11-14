function [x,y,z] = fitline_3dpts(x_sdata, y_sdata, z_sdata, axis_lims)

n_pts=100;

xyz = [x_sdata,y_sdata,z_sdata];
r0=mean(xyz);
xyz=bsxfun(@minus,xyz,r0);
[~,~,V]=svd(xyz,0);

t = linspace(axis_lims(1), axis_lims(2), n_pts);
x = r0(1) + t*V(1,1);
y = r0(2) + t*V(2,1);
z = r0(3) + t*V(3,1);

end