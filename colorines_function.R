##########################################################################################################################

###Function that overlaps the integrations position with the color domains to give a chromatin type to each of them###

##Input is the color domanis file and a data.frame containing at least the following columns:
    #Barcode
    #Chromosome
    #Position
#The color domains are load into the session and the file is modified to match the chromosome names

#Both are converted to a GRanges object and the match is done

##########################################################################################################################

colorines<- function(data.frame.name, color.domains){
    

    color.domains$seqname<- sub('chr','',color.domains$seqname)
    source("http://bioconductor.org/biocLite.R")
    library('GenomicRanges')

        color<- GRanges(seqnames= Rle(color.domains$seqname), ranges=IRanges(start= color.domains$start, end= color.domains$end), color= color.domains$chromatin)
        color<- narrow(color, start=3, end=-3)
        tmp<- GRanges(seqnames=Rle(data.frame.name$chr), ranges= IRanges(start= as.integer(data.frame.name$pos), end= as.integer(data.frame.name$pos)+1) )
        idx<- findOverlaps(tmp, color, select='first')
        integration.color<- color$color[idx]
    
    return(integration.color)
    }
