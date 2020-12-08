
clear all;
clc;
close all;

f_geno = '/home/aihual/meta_29112020/chr21.bim';

fileID = fopen(sprintf(f_geno));
bim_file = textscan(fileID,'%d %s %f %d %s %s');

fclose(fileID);

return
nsnp = length(bim_file{1});
snpid=bim_file{2};

snp_reduce=[3000 5000 8000];
snp_n=nsnp*ones(1,3)-snp_reduce;
snps=1:nsnp;
for i=1:size(snp_n,2)
    snp{i}=sort(randperm(nsnp,snp_n(i)));
    snp_w=fopen(['snp_',num2str(i),'.txt'],'w');
    
    for k=1:size(snp{i},2)
        fprintf(snp_w,'%s\n',snpid{snp{i}(k)});
    end
    fclose(snp_w);
end






