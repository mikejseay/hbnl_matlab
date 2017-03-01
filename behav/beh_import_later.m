function beh_table = beh_import( opt )
% collect behavioral data from a number of behavioral data .mats

% target the outpath
target_path = opt.outpath;

% index contents, find .m's
[mat_list, id_list] = coh_updatemats( target_path );

% indicate experiments to extract from each .mat
exps = opt.exps;

% indicate measures for each experiment
measures = {'acc', 'acc_withlate', 'medianrt', 'medianrt_withlate'};

% initialize data struct or w.e
% fields should be named by exp_descriptor_measure
beh_frame = struct('id_list', id_list );
for e = 1:length(exps)
    exp_ind = strcmpi( {opt.expstruct.name}, exps{e} );
    tmp_types = opt.expstruct(exp_ind).ttl_descriptors;
    for t = 1:length(tmp_types)
        for m = 1:length(measures)
            fname = [tmp_types{t}, '_', measures{m}];
            if ~isvarname(fname)
                if ~isempty( strfind(fname ,'$') )
                    fname( strfind(fname ,'$') ) = 'd';
                end
                if ~isempty( strfind(fname ,'+') )
                    fname( strfind(fname ,'+') ) = 'p';
                end
                if ~isempty( strfind(fname ,'-') )
                    fname( strfind(fname ,'-') ) = [];
                end
            end
            beh_frame.(fname) = [];
        end
    end
end

% actually import data into struct

for s = 1:length(mat_list);
    
    load( mat_list{s} )
    
    for e = 1:length(exps)
        if any( strcmpi( {output.name}, exps{e} ) )
            
            exp_ind = strcmpi( {output.name}, exps{e} );
            tmp_types = output(exp_ind).ttl_descriptors;
            
            for m = 1:length(measures)
                tmp_data = output(exp_ind).(measures{m});
                for t = 1:length(tmp_types)
                    fname = [tmp_types{t}, '_', measures{m}];
                    if ~isvarname(fname)
                        if ~isempty( strfind(fname ,'$') )
                            fname( strfind(fname ,'$') ) = 'd';
                        end
                        if ~isempty( strfind(fname ,'+') )
                            fname( strfind(fname ,'+') ) = 'p';
                        end
                        if ~isempty( strfind(fname ,'-') )
                            fname( strfind(fname ,'-') ) = [];
                        end
                    end
                    
                    beh_frame(s).(fname) = tmp_data(t);
                    
                end
            end
            
        end
    end
    
end

% remove fields that are completely nan-filled

beh_frame_fields = fieldnames(beh_frame);

for f = 1:length(beh_frame_fields)
    
    tmp_data = [ beh_frame.( beh_frame_fields{f} ) ];
    if all(isnan( tmp_data ))
        beh_frame = rmfield( beh_frame, beh_frame_fields{f});
    end
    
end

% convert to table and export to .csv

beh_table = struct2table(beh_frame);


end