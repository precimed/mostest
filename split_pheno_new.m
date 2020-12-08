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
Ncol = 10;
C_text = textscan(f1,formatSpec,Ncol);
format = repmat('%f ', [1 Ncol]);
pheno_file = textscan(f1,format);

fclose(fileID);
fclose(f1);

nsubj = length(fam_file{1});

fam_c1=fam_file{1};
fam_c2=fam_file{2};

subj_perm=randperm(nsubj);

%% distribute subj (fam file) into different sizes of cohorts
n1=2000;
n2=3000;
n3=nsubj-n1-n2;
cohort_num=[n1,n2,n3];
cohort_idex=[0,n1,n1+n2,n1+n2+n3];

n_cohort=size(cohort_num,2);
for i=1:n_cohort
    cohort{i}=sort(subj_perm(cohort_idex(i)+1:cohort_idex(i+1)));
    
    subj=fopen(['subj_',num2str(i),'.txt'],'w');
    pheno=fopen(['pheno_',num2str(i),'.txt'],'w');
    
    
    for k=1:cohort_num(i)
        fprintf(subj,'%s %s\n',fam_c1{cohort{i}(k)},fam_c2{cohort{i}(k)});
        
        for j=1:10
            fprintf(pheno,'%f ',pheno_file{j}(cohort{i}(k)));
            if j==10
                fprintf(pheno,'\n');
            end
        end
    end
    fclose(subj);
    fclose(pheno);
end
%intersect(cohort{1},cohort{2}) %check no overlapping



