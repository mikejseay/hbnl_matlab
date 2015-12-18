%% determine dependencies

function_list={'coh_createopt.m','coh_makesuperbatches.m','list2cell.m',...
    'basename.m','coh_calc.m','full_load.m','pop_jointprob.m','jointprob.m'};

deps=matlab.codetools.requiredFilesAndProducts(function_list)';