clear all;
clc;
close all;

f_subject = '/home/aihual/meta_02112020/chr21.fam'

fileID = fopen(sprintf(f_subject));
fam_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);

nsubj = length(fam_file{1});


subj_c1=fam_file{1};
subj_c2=fam_file{2};


subj_perm=randperm(nsubj);

n_cohort=4;

subj_sub=reshape(subj_perm,[nsubj/n_cohort,n_cohort]);

cohort_sub=sort(subj_sub);



cohort_1=sort(subj_sub(:,1));

cohort_2=sort(subj_sub(:,2));



%%C=intersect(cohort_1, cohort_2 )%% double check if cohort_1 and cohort_2
%%are non-overlapping

subj_1 = fopen('cohort_1.txt','w');

subj_2 = fopen('cohort_2.txt','w');

for k=1:size(cohort_sub,1)
    fprintf(subj_1,'%s %s\n',subj_c1{cohort_1(k)},subj_c2{cohort_1(k)});
    fprintf(subj_2,'%s %s\n',subj_c1{cohort_2(k)},subj_c2{cohort_2(k)});
    
end

fclose(subj_1);
fclose(subj_2);


f_pheno = '/home/aihual/mostest_02112020/pheno.txt'

f1 = fopen(sprintf(f_pheno));
return
s_file = textscan(f1,'%s %s %s %s %s %s %s %s %s %s');
fclose(f1);


pheno_1 = fopen('subj_1.txt','w');

pheno_2 = fopen('subj_2.txt','w');

for k=1:size(cohort_1,1)
    for j=1:10
        fprintf(pheno_1,'%s\t',s_file{j}{cohort_1(k)+1});
        if j==10
            fprintf(pheno_1,'\n');
        end
    end
    
end



fclose(pheno_1);
fclose(pheno_2);



