
clear all;
clc;
close all;

f_geno = '/home/aihual/meta_29112020/chr21.bim';

fileID = fopen(sprintf(f_geno));
bim_file = textscan(fileID,'%d %s %f %d %s %s');

fclose(fileID);

nsnp = length(bim_file{1});
snpid=bim_file{2};

max_reduce=floor(nsnp/60);% set a up limit for the number of dropped snps of each cohort
snp_reduce=[randi(max_reduce,1) randi(max_reduce,1) randi(max_reduce,1)];

num_cohort=3;
snp_n=nsnp*ones(1,num_cohort)-snp_reduce;% snp list for each cohort
snps=1:nsnp;
for i=1:num_cohort
    snp{i}=sort(randperm(nsnp,snp_n(i)));
    snp_w=fopen(['snp_',num2str(i),'.txt'],'w');
    
    for k=1:snp_n(i)
        fprintf(snp_w,'%s\n',snpid{snp{i}(k)});
    end
    fclose(snp_w);
end






