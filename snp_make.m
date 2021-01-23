clear all;
close all;
clc;


f1=fopen('chr21.bim','r');
s=textscan(f1,'%d%s%f%d%s%s\n');
fclose(f1);
s1=[s{1}];%no.chr
s2=s{2};%snp names


s1(1:10)=0;

nc=1000;%%no.variants

[idx,val]=find(s1==21);%find index of chr21 

id1=idx(1);
id2=idx(end);
no_chr21=id2-id1+1;%no.chr21
snp_rand=randperm(no_chr21,nc);%randomly take nc=1000 unrepeated numbers from [1:no_chr21]
snp_eff=id1+sort(snp_rand)-1;%index of snps with nonzeros
save('snp_eff.mat','snp_eff')


beta_mu=0;
beta_sigma=0.3;

no_traits=10;

for i=1:no_traits
    s3=zeros(size(s1,1),1);%vector of beta

    % beta, normal distribution
    s3(snp_eff)=normrnd(beta_mu,beta_sigma,nc,1);
    
    fid(i) = fopen(['snp_',num2str(i),'.txt'],'w');
  
    for k=1:size(snp_eff,2)
        %fprintf(fid(i),'%s\t%f\r\n',s2{k,:},s3(k,:));
       fprintf(fid(i),'%s\t%f\n',s2{snp_eff(k),:},s3(snp_eff(k),:));
    end
    fclose(fid(i));
end
