function genomat = PlinkRead_binary2(nsubj,snps,fileprefix)
persistent geno_values

% Written by Chun 2015
if ~issorted(snps), error('PlinkRead_binary2 expect sorted list of snps'); end;
nsnp = length(snps);

% bit shift to generate the genovalue matrix
bedprefix = sprintf('%s.bed', fileprefix);

%load('PlinkRead_binary2_geno_values.mat')
if isempty(geno_values)
    geno_values = [
     2    2    2    2
    -1    2    2    2
     1    2    2    2
     0    2    2    2
     2   -1    2    2
    -1   -1    2    2
     1   -1    2    2
     0   -1    2    2
     2    1    2    2
    -1    1    2    2
     1    1    2    2
     0    1    2    2
     2    0    2    2
    -1    0    2    2
     1    0    2    2
     0    0    2    2
     2    2   -1    2
    -1    2   -1    2
     1    2   -1    2
     0    2   -1    2
     2   -1   -1    2
    -1   -1   -1    2
     1   -1   -1    2
     0   -1   -1    2
     2    1   -1    2
    -1    1   -1    2
     1    1   -1    2
     0    1   -1    2
     2    0   -1    2
    -1    0   -1    2
     1    0   -1    2
     0    0   -1    2
     2    2    1    2
    -1    2    1    2
     1    2    1    2
     0    2    1    2
     2   -1    1    2
    -1   -1    1    2
     1   -1    1    2
     0   -1    1    2
     2    1    1    2
    -1    1    1    2
     1    1    1    2
     0    1    1    2
     2    0    1    2
    -1    0    1    2
     1    0    1    2
     0    0    1    2
     2    2    0    2
    -1    2    0    2
     1    2    0    2
     0    2    0    2
     2   -1    0    2
    -1   -1    0    2
     1   -1    0    2
     0   -1    0    2
     2    1    0    2
    -1    1    0    2
     1    1    0    2
     0    1    0    2
     2    0    0    2
    -1    0    0    2
     1    0    0    2
     0    0    0    2
     2    2    2   -1
    -1    2    2   -1
     1    2    2   -1
     0    2    2   -1
     2   -1    2   -1
    -1   -1    2   -1
     1   -1    2   -1
     0   -1    2   -1
     2    1    2   -1
    -1    1    2   -1
     1    1    2   -1
     0    1    2   -1
     2    0    2   -1
    -1    0    2   -1
     1    0    2   -1
     0    0    2   -1
     2    2   -1   -1
    -1    2   -1   -1
     1    2   -1   -1
     0    2   -1   -1
     2   -1   -1   -1
    -1   -1   -1   -1
     1   -1   -1   -1
     0   -1   -1   -1
     2    1   -1   -1
    -1    1   -1   -1
     1    1   -1   -1
     0    1   -1   -1
     2    0   -1   -1
    -1    0   -1   -1
     1    0   -1   -1
     0    0   -1   -1
     2    2    1   -1
    -1    2    1   -1
     1    2    1   -1
     0    2    1   -1
     2   -1    1   -1
    -1   -1    1   -1
     1   -1    1   -1
     0   -1    1   -1
     2    1    1   -1
    -1    1    1   -1
     1    1    1   -1
     0    1    1   -1
     2    0    1   -1
    -1    0    1   -1
     1    0    1   -1
     0    0    1   -1
     2    2    0   -1
    -1    2    0   -1
     1    2    0   -1
     0    2    0   -1
     2   -1    0   -1
    -1   -1    0   -1
     1   -1    0   -1
     0   -1    0   -1
     2    1    0   -1
    -1    1    0   -1
     1    1    0   -1
     0    1    0   -1
     2    0    0   -1
    -1    0    0   -1
     1    0    0   -1
     0    0    0   -1
     2    2    2    1
    -1    2    2    1
     1    2    2    1
     0    2    2    1
     2   -1    2    1
    -1   -1    2    1
     1   -1    2    1
     0   -1    2    1
     2    1    2    1
    -1    1    2    1
     1    1    2    1
     0    1    2    1
     2    0    2    1
    -1    0    2    1
     1    0    2    1
     0    0    2    1
     2    2   -1    1
    -1    2   -1    1
     1    2   -1    1
     0    2   -1    1
     2   -1   -1    1
    -1   -1   -1    1
     1   -1   -1    1
     0   -1   -1    1
     2    1   -1    1
    -1    1   -1    1
     1    1   -1    1
     0    1   -1    1
     2    0   -1    1
    -1    0   -1    1
     1    0   -1    1
     0    0   -1    1
     2    2    1    1
    -1    2    1    1
     1    2    1    1
     0    2    1    1
     2   -1    1    1
    -1   -1    1    1
     1   -1    1    1
     0   -1    1    1
     2    1    1    1
    -1    1    1    1
     1    1    1    1
     0    1    1    1
     2    0    1    1
    -1    0    1    1
     1    0    1    1
     0    0    1    1
     2    2    0    1
    -1    2    0    1
     1    2    0    1
     0    2    0    1
     2   -1    0    1
    -1   -1    0    1
     1   -1    0    1
     0   -1    0    1
     2    1    0    1
    -1    1    0    1
     1    1    0    1
     0    1    0    1
     2    0    0    1
    -1    0    0    1
     1    0    0    1
     0    0    0    1
     2    2    2    0
    -1    2    2    0
     1    2    2    0
     0    2    2    0
     2   -1    2    0
    -1   -1    2    0
     1   -1    2    0
     0   -1    2    0
     2    1    2    0
    -1    1    2    0
     1    1    2    0
     0    1    2    0
     2    0    2    0
    -1    0    2    0
     1    0    2    0
     0    0    2    0
     2    2   -1    0
    -1    2   -1    0
     1    2   -1    0
     0    2   -1    0
     2   -1   -1    0
    -1   -1   -1    0
     1   -1   -1    0
     0   -1   -1    0
     2    1   -1    0
    -1    1   -1    0
     1    1   -1    0
     0    1   -1    0
     2    0   -1    0
    -1    0   -1    0
     1    0   -1    0
     0    0   -1    0
     2    2    1    0
    -1    2    1    0
     1    2    1    0
     0    2    1    0
     2   -1    1    0
    -1   -1    1    0
     1   -1    1    0
     0   -1    1    0
     2    1    1    0
    -1    1    1    0
     1    1    1    0
     0    1    1    0
     2    0    1    0
    -1    0    1    0
     1    0    1    0
     0    0    1    0
     2    2    0    0
    -1    2    0    0
     1    2    0    0
     0    2    0    0
     2   -1    0    0
    -1   -1    0    0
     1   -1    0    0
     0   -1    0    0
     2    1    0    0
    -1    1    0    0
     1    1    0    0
     0    1    0    0
     2    0    0    0
    -1    0    0    0
     1    0    0    0
     0    0    0    0
    ];
    geno_values = cast(geno_values, 'int8');
end

if isempty(geno_values)
    geno_values = zeros(256,4,'int8');
    geno_code = [-1,1,0,2];
    shiftind = [0,2,4,6];
    indvec=zeros(1,4);

    for i = 1:256;
        ishift = int16(i-1);
        for j = 1:4;
            indvec(j) = bitand(bitsra(ishift,shiftind(j)),3) ;
        end
        indvec(indvec == 0) = 4;
        geno_values(i,:) = geno_code(indvec);
    end
    %save('PlinkRead_binary2_geno_values.mat', '-v7', 'geno_values')
end

% Read in the binary file
bedid = fopen(bedprefix);
genobin = uint16(fread(bedid, 3));

% Check magic number
if genobin(1) ~= 108;
	error('- Not a valid Plink BED file \r\n');
elseif genobin(2) ~= 27;
	error('- Not a valid Plink BED file \r\n');
elseif genobin(3) ~= 1;
	error('- Not in SNP-major format \r\n');
end

n_bytes = ceil(nsubj/4);
genomat = zeros(nsubj,nsnp,'int8');
for i = 1:nsnp;
    fseek(bedid, 3 + (snps(i) - 1) * n_bytes, 'bof');
    genobin = uint16(fread(bedid, n_bytes));
    if length(genobin) ~= n_bytes, error('-- Invalid number of entries from the bed \r\n'); end
    tmp_values = geno_values(genobin + 1, :)';
    tmp_values = tmp_values(:);
    genomat(:,i) = tmp_values(1:nsubj);
end
fclose(bedid);

end