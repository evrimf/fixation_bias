##function to remove lines that violate HWE

import math
import scipy
from scipy.stats import chi2

def HW(population_info, infile, outfile):

    population_info=open(population_info,"r")
    inf=open(infile,"r")
    outf=open(outfile,"a")

    sub_or_super=2 # 2 for using superpopulations, 1 for subpopulations as sample grouping

    number_of _columns_with_bed_info=5
    number_of_populations_allow_to_violate=1
    pval_treshold=0.01

    pop_info={}
    count=0
    population_info.readline()  #skip first line

    for line in population_info:
        line=line.split()
        try:
            pop_info[line[sub_or_super]] += [count] # get indices as (list) value for each population code as key

        except KeyError:
            pop_info[line[sub_or_super]]=[count]

    population_info.close()
    
# within each line, for each population
# count genotypes, calculate allele frequencies and hw expected genotypes
# calculate chi square and pval
# if pval is smaller than threshold (violates hwe) for given # of populations
# do nothing and continue with next line
# else, print bed info columns and global reference allele frequency

    for line in inf:
        violations = 0
        line = line.split()
        bed_info = line[0:number_of_lines_with_bed_info]
        line = line[number_of_lines_with_bed_info:]
        for population in pop_info.keys():
            pop_columns = [line[i] for i in pop_info[population]]  # select genotypes belonging to current population

            count_refhomo = float(pop_columns.count("0|0"))
            count_hetero = float(pop_columns.count("0|1") + pop_columns.count("1|0"))
            count_althomo = float(pop_columns.count("1|1"))
            N = float(len(pop_columns))

            allele_freq_ref = count_refhomo/N + count_hetero/N/2
            allele_freq_alt = count_althomo/N + count_hetero/N/2

            if (allele_freq_ref == 0) | (allele_freq_ref == 1):
                continue # these can not violate hwe

            exp_refhomo = N * allele_freq_ref**2
            exp_hetero = 2 * N * allele_freq_ref * allele_freq_alt
            exp_althomo = N * allele_freq_alt**2

            chisqx = (count_refhomo - exp_refhomo)**2/exp_refhomo + (count_hetero - exp_hetero)**2/exp_hetero + (count_althomo - exp_althomo)**2/exp_althomo

            pvalx = 1 - chi2.cdf(chisqx, 1)

            if pvalx < pval_threshold:
                violations += 1

        if violations>number_of_populations_allow_to_violate:
            continue

		else:
		
			count_refhomo_g = float(line.count("0|0"))                   #0 referans
			count_hetero_g = float(line.count("0|1") + line.count("1|0"))
			Numberof_genotypes = float(len(line))
		
			allele_freq_ref_g = count_refhomo_g/N_g + count_hetero_g/N_g/2

#		print "count_refhomo_g", count_refhomo_g
#		print "count_hetero_g", count_hetero_g
#		print "N_g", N_g
#		print "allele_freq_ref_g", allele_freq_ref_g

			output_line = '\t'.join(bed_info + [str(allele_freq_ref_g)])+'\n'
			outf.write(output_line)

	inf.close()
    outf.close()
         
    
