xpost_2_atad = load('posteriorModel_2_ATAD.mat');
x_post_2_atad = xpost_2_atad.results{end,1};
mv_post_2_atad = mean(10.^x_post_2_atad(10000:end,end-1:end));
vv_post_2_atad = var(10.^x_post_2_atad(10000:end,end-1:end));

xpost_2_rad = load('posteriorModel_2_RAD.mat');
x_post_2_rad = xpost_2_rad.results{end,1};
mv_post_2_rad = mean(10.^x_post_2_rad(10000:end,end-1:end));
vv_post_2_rad = var(10.^x_post_2_rad(10000:end,end-1:end));