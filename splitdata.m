clear all;
clc;
close all;


%% randomly generate n_cohort geno data
f_geno = '/home/aihual/meta_29112020/chr21.bim';

fileID = fopen(sprintf(f_geno));
bim_file = textscan(fileID,'%d %s %f %d %s %s');

fclose(fileID);

nsnp = length(bim_file{1});
snpid=bim_file{2};
n_cohort=3;
max_reduce=floor(nsnp/100);% set a up limit for the number of dropped snps of each cohort
%snp_reduce=[randi(max_reduce,1) randi(max_reduce,1) randi(max_reduce,1)];
for n=1:n_cohort
    snp_reduce(n)=randi(max_reduce,1);
end

snp_n=nsnp*ones(1,n_cohort)-snp_reduce;% snp list for each cohort
snps=1:nsnp;
for i=1:n_cohort
    snp{i}=sort(randperm(nsnp,snp_n(i)));
    snp_w=fopen(['snps_',num2str(i),'.txt'],'w');
    
    for k=1:snp_n(i)
        fprintf(snp_w,'%s\n',snpid{snp{i}(k)});
    end
    fclose(snp_w);
end

%% divide .fam and pheno.txt into n_cohort
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
%[cohort{1},cohort{2},cohort{3}] = dividerand(nsubj,0.4,0.3,0.3);
cohort=divide(nsubj,n_cohort);

for i=1:n_cohort
    subj=fopen(['subj_',num2str(i),'.txt'],'w');
    pheno=fopen(['cohort_',num2str(i),'.txt'],'w');
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


