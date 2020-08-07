get_bed_geno = function(bfile, isnp, nsamples, byte_geno_map=NULL) {
    # Args:
    #   bfile:    plink bim/fam/bed file prefix
    #   isnp:     vector of SNPs indices to extract
    #   nsampels: number of samples in bfile
    #   byte_geno_map: 2D array of byte to genotypes mapping, produced with
    #                  get_byte_geno_map function. This argument is optional,
    #                  but providing it will improve runtime if the function
    #                  is called multiple times.
    # Return:
    #   genomat: [nsamples x nspns] matrix of genotypes
    #
    # Call example:
    #   Extract genotypes for first 100 SNPs from '/path/to/chr1.bed' file
    #   with 10000 individuals:
    #   genomat = get_bed_geno('/path/to/chr1',1:100,10000)

    if( is.null(byte_geno_map) ) {
        byte_geno_map = get_byte_geno_map()
    }

    bed_path = paste(bfile, 'bed', sep='.')
    bedf = file(bed_path, 'rb')

    # check magic bytes
    magic_bytes = readBin(bedf, 'integer', n=3, size=1, signed=FALSE)
    if( !all( magic_bytes == c(108, 27, 1)) ) {
        stop('Not a valid Plink BED file')
    }

    genomat = matrix(0, nrow=nsamples, ncol=length(isnp))
    nbytes = if( nsamples%%4 > 0 ) nsamples%/%4 + 1 else nsamples%/%4
    j = 1
    for( i in isnp ) {
        seek(bedf, 3 + (i-1)*nbytes)
        bvec = readBin(bedf, 'integer', n=nbytes, size=1, signed=FALSE)
        geno_i = byte_geno_map[,bvec + 1]
        geno_i = matrix(geno_i, nrow=4*nbytes, ncol=1) # reshape to column vector
        geno_i = geno_i[1:nsamples]
        genomat[,j] = geno_i
        j = j + 1
    }

    close(bedf)

    return(genomat)
}


get_byte_geno_map = function() {
    # Return:
    #   byte_geno_map: [4 x 256] array, byte_geno_map[,i] = [a,b,c,d], where
    #                  a,b,c,d are from {0,1,2,-1}, genotypes encodede by byte i.
    #                  See definition of Plink's bed format: https://www.cog-genomics.org/plink2/formats#bed
    geno_codes = c(2, -1, 1, 0)
    bit_shifts = c(0, 2, 4, 6)
    byte_geno_map = matrix(0, nrow=4, ncol=256)
    for( b in 0:255 ) {
        byte_geno_map[,b+1] = geno_codes[bitwAnd(bitwShiftR(b, bit_shifts), 3) + 1]
    }
    return(byte_geno_map)
}

