files2merge = ["/mnt/seagate10/projects/mostest_vertexwise/results/experiment8/lh_rh.pial.0.fsaverage3.reg_10.stats.mat", ...
               "/mnt/seagate10/projects/mostest_vertexwise/results/experiment8/lh_rh.pial.0.fsaverage3.reg_10.stats_1.mat"];

out_path = "/mnt/seagate10/projects/mostest_vertexwise/results/experiment8/lh_rh.pial.0.fsaverage3.reg_10.stats.stats_1.mat";

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

