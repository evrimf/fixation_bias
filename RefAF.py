#!/usr/bin/python
###function to find reference allele frequencies after HW filtering

def RefAF(infile,outfile):
    
    inf=open(infile, "r")
    outf=open(outfile,"a")

    number_of _columns_with_bed_info=5

    for line in inf:

        line=line.split()
        bed_info=line[0:number_of_columns_with_bed_info]
        line=line[number_of_columns_with_bed_info:]

        count_refhomo_g=float(line.count("0|0"))
        count_hetero_g=float(line.count("0|1")+line.count("1|0"))
        N_g=float(len(line))

        allele_freq_ref_g=count_refhomo_g/N_g + count_hetero_g/N_g/2

        output_line='\t'.join(bed_info+[str(allele_freq_ref_g)])+'\n'
        outf.write(output_line)

    inf.close()
    outf.close()

    
    
    
    
