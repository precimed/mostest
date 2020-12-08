clear all;
clc;
close all;

%% fetch fam file and pheno file
f_subject = '/home/aihual/meta_29112020/chr21.fam';
fileID = fopen(sprintf(f_subject));
fam_file = textscan(fileID,'%s %s %s %s %s %s');

f_pheno = '/home/aihual/meta_29112020/pheno.txt';
f1 = fopen(sprintf(f_pheno));
formatSpec = '%s';
Ncol = 10;% no.columns of pheno
C_text = textscan(f1,formatSpec,Ncol);
format = repmat('%f ', [1 Ncol]);
pheno_file = textscan(f1,format);

fclose(fileID);
fclose(f1);

nsubj = length(fam_file{1});

fam_c1=fam_file{1};
fam_c2=fam_file{2};

%% distribute subj (fam file) into different sizes of cohorts
n_cohort=3;
[cohort{1},cohort{2},cohort{3}] = dividerand(nsubj,0.4,0.3,0.3);

for i=1:n_cohort
    subj=fopen(['subj_',num2str(i),'.txt'],'w');
    pheno=fopen(['pheno_',num2str(i),'.txt'],'w');
    cohort_num=size(cohort{i},2);
    
    for k=1:cohort_num
        fprintf(subj,'%s %s\n',fam_c1{cohort{i}(k)},fam_c2{cohort{i}(k)});
        
        for j=1:Ncol
            fprintf(pheno,'%f ',pheno_file{j}(cohort{i}(k)));
            if j==Ncol
                fprintf(pheno,'\n');
            end
        end
    end
    fclose(subj);
    fclose(pheno);
end
%intersect(cohort{1},cohort{2}) %check no overlapping



