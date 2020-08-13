files2merge = ["/mnt/seagate10/projects/mostest_vertexwise/results/experiment_roi/stats/thickness/lh_rh.thickness.UKBB.no_avg.reg_8.stats_mmil1.mat", ...
"/mnt/seagate10/projects/mostest_vertexwise/results/experiment_roi/stats/thickness/lh_rh.thickness.UKBB.no_avg.reg_8.stats_mmil2.mat", ...
"/mnt/seagate10/projects/mostest_vertexwise/results/experiment_roi/stats/thickness/lh_rh.thickness.UKBB.no_avg.reg_8.stats_mmil3.mat", ...
"/mnt/seagate10/projects/mostest_vertexwise/results/experiment_roi/stats/thickness/lh_rh.thickness.UKBB.no_avg.reg_8.stats_mmil4.mat", ...
"/mnt/seagate10/projects/mostest_vertexwise/results/experiment_roi/stats/thickness/lh_rh.thickness.UKBB.no_avg.reg_8.stats_mmil5.mat", ...
"/mnt/seagate10/projects/mostest_vertexwise/results/experiment_roi/stats/thickness/lh_rh.thickness.UKBB.no_avg.reg_8.stats_mmil6.mat"];

out_path = "/mnt/seagate10/projects/mostest_vertexwise/results/experiment_roi/stats/thickness/lh_rh.thickness.UKBB.no_avg.reg_8.no_avg.stats_mmil1.stats_mmil2.stats_mmil3.stats_mmil4.stats_mmil5.stats_mmil6.mat";

fprintf('Loading %s ... ', files2merge(1));
load(files2merge(1));
fprintf('Done\n');
merge_mostvecs_perm = mostvecs_perm;
merge_minpvecs_perm = minpvecs_perm;

for i = 2:numel(files2merge)
    fprintf('Loading %s ... ', files2merge(i));
    load(files2merge(i));
    fprintf('Done\n');
    merge_mostvecs_perm = [merge_mostvecs_perm, mostvecs_perm];
    merge_minpvecs_perm = [merge_minpvecs_perm, minpvecs_perm];
end

mostvecs_perm = merge_mostvecs_perm;
minpvecs_perm = merge_minpvecs_perm;

[nsnp nperm] = size(minpvecs_perm);
fprintf('%d permutations with %d SNPs\n', nperm, nsnp);

fprintf('saving %s as -v7.3 ... ', out_path);
save(out_path, '-v7.3', 'nvec', 'freqvec', 'gwas_time_sec', 'mostvecs_orig', 'minpvecs_orig', 'mostvecs_perm', 'minpvecs_perm');
fprintf('Done\n')

