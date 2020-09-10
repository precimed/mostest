function geno_mat = bed2mat(bfile, snp_idx)
    % Args:
    %   bfile:   plink bfile, prefix of bed/bim/fam files.
    %   snp_idx: array of SNP indices to fetch from bfile.bed.
    %
    % Return:
    %   geno_mat: nsampels x numel(snp_idx) matrix, nsampels = number of samples in bfile.fam.

    persistent byte_geno_map;
    if isempty(byte_geno_map)
        geno_codes = [2,-1,1,0];
        bit_shifts = [0,-2,-4,-6];
        byte_geno_map = zeros(4,256,'int8');
        for i = 0:255
            byte_geno_map(:,i+1) = geno_codes(bitand(bitshift(i, bit_shifts), 3) + 1);
        end
    end

    persistent bed last_bfile nsnps nsamples;
    bed_path = sprintf('%s.bed', bfile);

    if isempty(last_bfile) || ~strcmp(last_bfile, bfile)
        last_bfile = bfile;
        % check magic bits
        bed_id = fopen(bed_path);
        a = uint8(fread(bed_id, 3));
        if ~all(a == [108;27;1])
            error('%s is not a valid plink bed file\n', bed_path);
        end
        fclose(bed_id);

        bim_id = fopen(sprintf('%s.bim', bfile));
        bim_file = textscan(bim_id,'%s %*[^\n]');
        fclose(bim_id);
        nsnps = length(bim_file{1});

        fam_id = fopen(sprintf('%s.fam', bfile));
        fam_file = textscan(fam_id,'%s %*[^\n]');
        fclose(fam_id);
        nsamples = length(fam_file{1});

        nrows = ceil(nsamples/4);
        fprintf('loading %s\n', bed_path);
        bed = memmapfile(bed_path,'Offset',3,'Format',{'uint8',[nrows,nsnps],'geno'});

        fprintf('%d snps, %d samples loaded from %s\n', nsnps, nsamples, bfile);
    end

    if max(snp_idx) > nsnps
        error('Max queried SNP index is > than the number of SNPs in %s\n', bfile);
    end
   
    ui16 = uint16(bed.Data.geno(:,snp_idx)) + 1;
    geno_mat = reshape(byte_geno_map(:,ui16), 4*nrows, numel(snp_idx));
    geno_mat = geno_mat(1:nsamples,:);
end
