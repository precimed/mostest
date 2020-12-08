
f_geno = '/home/aihua/mostest_debug/chr21.bim'

fileID = fopen(sprintf(f_geno));
bim_file = textscan(fileID,'%d %s %f %d %s %s');
fclose(fileID);

nsnp = length(bim_file{1});
snpid=bim_file{2};

unique_1=30;
unique_2=40;
overlap_n=nsnp-unique_1-unique_2-50;
snp_t=1:nsnp;


snp_unique_1=randperm(nsnp,unique_1);
snp_left_1=setdiff(snp_t,snp_unique_1);

index_unique_2=randperm(nsnp-unique_1, unique_2);
snp_unique_2=snp_left_1(index_unique_2);

snp_left_2=setdiff(snp_left_1,snp_unique_2);

index_overlap=randperm(nsnp-unique_1-unique_2,overlap_n);
snp_overlap=snp_left_2(index_overlap);

snp_1=sort(horzcat(snp_unique_1,snp_overlap));
snp_2=sort(horzcat(snp_unique_2,snp_overlap));

snpid_1 = fopen('snpid_1.txt','w');

snpid_2 = fopen('snpid_2.txt','w');

for k=1:size(snp_1,1)
    fprintf(snpid_1,'%s\n',snpid{snp_1});
end

for k=1:size(snp_2,1)
    fprintf(snpid_2,'%s\n',snpid{snp_2});
end


fclose(snpid_1);
fclose(snpid_2);




