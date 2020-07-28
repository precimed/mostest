stepwise_features_dir = "/mnt/seagate10/projects/mostest_vertexwise/pheno/stepwise_features";
ids_file = sprintf('%s/vertex_subjectlist', stepwise_features_dir);
covar_file = sprintf('%s/demographics.mostest.csv', stepwise_features_dir);

fwhm_arr = [0,20];
feature_arr = ["pial","thickness"];
hemi_arr = ["lh_rh"];

for i_fwhm=1:numel(fwhm_arr)
    fwhm = fwhm_arr(i_fwhm);
    for i_feature=1:numel(feature_arr)
        feature = feature_arr(i_feature);
        for i_hemi=1:numel(hemi_arr)
            hemi = hemi_arr(i_hemi);
            pheno_file = sprintf('%s/fwhm%d/%s.%s.%d.fsaverage3.csv', stepwise_features_dir, fwhm, hemi, feature, fwhm);
            fname_out = sprintf('%s/fwhm%d/%s.%s.%d.fsaverage3.resid.mat', stepwise_features_dir, fwhm, hemi, feature, fwhm);
            if exist(pheno_file, 'file') == 2
                fprintf('%s\n', pheno_file);
                pheno_mat = resid(pheno_file, ids_file, covar_file);
                save(fname_out, '-v7.3', 'pheno_mat');
                fprintf('    %s saved\n', fname_out);
            else
                fprintf('%s not found\n', pheno_file);
            end
        end
    end
end
