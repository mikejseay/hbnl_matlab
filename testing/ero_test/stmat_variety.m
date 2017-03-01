v4_mc16_21 = '/processed_data/mat-files-v40/mc16-21/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_1_a1_40001002_256_cnt.h1.a.st.mat';
v4_mc16_32 = '/processed_data/mat-files-v40/mc16-32/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_1_a1_40001007_256_cnt.h1.a.st.mat';
v4_mc16_64 = '/processed_data/mat-files-v40/mc16-64/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_2_a1_40001013_256_cnt.h1.a.st.mat';
v4_ns16_32 = '/processed_data/mat-files-v40/ns16-32/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_3_a1_40188024_256_cnt.h1.a.st.mat';
v4_ns16_64 = '/processed_data/mat-files-v40/ns16-64/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_3_a1_40025004_256_cnt.h1.a.st.mat';
v4_ns32_64 = '/processed_data/mat-files-v40/ns32-64/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_5_a1_40001015_256_cnt.h1.a.st.mat';

v6_mc16_21 = '/processed_data/mat-files-v60/mc16-21/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_1_a1_40001002_256_cnt.h1.a.st.mat';
v6_mc16_32 = '/processed_data/mat-files-v60/mc16-32/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_1_a1_40001007_256_cnt.h1.a.st.mat';
v6_mc16_64 = '/processed_data/mat-files-v60/mc16-64/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_2_a1_40001013_256_cnt.h1.a.st.mat';
v6_ns16_32 = '/processed_data/mat-files-v60/ns16-32/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_3_a1_40188024_256_cnt.h1.a.st.mat';
v6_ns16_64 = '/processed_data/mat-files-v60/ns16-64/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_3_a1_40025004_256_cnt.h1.a.st.mat';
v6_ns32_64 = '/processed_data/mat-files-v60/ns32-64/ant/ant-a/e1-n10-s9-t100-v800/suny/ant_5_a1_40001015_256_cnt.h1.a.st.mat';

example_stmats = {v4_mc16_21, v4_mc16_32, v4_mc16_64, v4_ns16_32, v4_ns16_64, ...
    v4_ns32_64, v6_mc16_21, v6_mc16_32, v6_mc16_64, v6_ns16_32, v6_ns16_64, ...
    v6_ns32_64};

for ind=1:12
    disp(example_stmats{ind})
    load(example_stmats{ind})
    check = 0;
end