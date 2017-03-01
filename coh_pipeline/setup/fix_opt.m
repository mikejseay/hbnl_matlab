function opt=fix_opt(param_struct,outpath,coordsfile)

opt=param_struct;
opt.wavelet_scales = opt.scale;
opt = rmfield(opt, 'scale');
opt.outpath=outpath;
opt.coords_file=coordsfile;


end