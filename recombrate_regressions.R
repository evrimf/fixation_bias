##2017-07-25
##Logistic regressions with DECODE recombination rates
#decode_table -> DECODE recombination rates
decode_table=read.table("/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/processed_data/pattern_recomb_fixed-poly_genic-nongenic_score.bed")
colnames(decode_table)=c("chr","start pos.","stop pos.","mut.pattern","recombination rate","fixed|poly","genic|nongenic","hscore")
decode_table=decode_table[,4:8]
head(decode_table)

##all mutation types
##factors separate and with interaction
decode_regression_recomb=glm(decode_table[,"fixed|poly"]~decode_table[,"recombination rate"],family="binomial")
decode_regression_genic=glm(decode_table[,"fixed|poly"]~factor(decode_table[,"genic|nongenic"]),family="binomial")
decode_regression_hscore=glm(decode_table[,"fixed|poly"]~decode_table[,"hscore"],family="binomial")
decode_regression_recomb_hscore=glm(decode_table[,"fixed|poly"]~decode_table[,"recombination rate"]*decode_table[,"hscore"],family="binomial")
decode_regression_all=glm(decode_table[,"fixed|poly"]~decode_table[,"recombination rate"]+factor(decode_table[,"genic|nongenic"])+decode_table[,"hscore"]+decode_table[,"recombination rate"]*decode_table[,"hscore"],family="binomial")

summary(decode_regression_recomb)
summary(decode_regression_genic)
summary(decode_regression_hscore)
summary(decode_regression_recomb_hscore)
summary(decode_regression_all)

save(decode_regression_recomb,decode_regression_genic,decode_regression_hscore,decode_regression_recomb_hscore,decode_regression_all,file="/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/rdata/decode_recomb_regression_all.Rdata")

##betwen C<->G
cg_indices=grep("[A-Z]{2}[C|G][A-Z]{2}->[A-Z]{2}[C|G][A-Z]{2}",decode_table[,1])

decode_cg_regression_recomb=glm(decode_table[cg_indices,"fixed|poly"]~decode_table[cg_indices,"recombination rate"],family="binomial")
decode_cg_regression_genic=glm(decode_table[cg_indices,"fixed|poly"]~factor(decode_table[cg_indices,"genic|nongenic"]),family="binomial")
decode_cg_regression_hscore=glm(decode_table[cg_indices,"fixed|poly"]~decode_table[cg_indices,"hscore"],family="binomial")
decode_cg_regression_recomb_hscore=glm(decode_table[cg_indices,"fixed|poly"]~decode_table[cg_indices,"recombination rate"]*decode_table[cg_indices,"hscore"],family="binomial")
decode_cg_regression_all=glm(decode_table[cg_indices,"fixed|poly"]~decode_table[cg_indices,"recombination rate"]+factor(decode_table[cg_indices,"genic|nongenic"])+decode_table[cg_indices,"hscore"]+decode_table[cg_indices,"recombination rate"]*decode_table[cg_indices,"hscore"],family="binomial")

summary(decode_cg_regression_recomb)
summary(decode_cg_regression_genic)
summary(decode_cg_regression_hscore)
summary(decode_cg_regression_recomb_hscore)
summary(decode_cg_regression_all)

save(decode_cg_regression_recomb,decode_cg_regression_genic,decode_cg_regression_hscore,decode_cg_regression_recomb_hscore,decode_cg_regression_all,file="/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/rdata/decode_recomb_regression_cg.Rdata")

##between AT->CG patterns (gBGC control)
##output:/Users/evrimfer/Desktop/Mutation_Fixation_Bias/results/regression_results/2017-07-25-decode_gbgc_regression.txt

gbgc_indices=grep("[A-Z]{2}[A|T][A-Z]{2}->[A-Z]{2}[C|G][A-Z]{2}", decode_table[,1])

decode_gbgc_regression_recomb=glm(decode_table[gbgc_indices,"fixed|poly"]~decode_table[gbgc_indices,"recombination rate"],family="binomial")
decode_gbgc_regression_genic=glm(decode_table[gbgc_indices,"fixed|poly"]~factor(decode_table[gbgc_indices,"genic|nongenic"]),family="binomial")
decode_gbgc_regression_hscore=glm(decode_table[gbgc_indices,"fixed|poly"]~decode_table[gbgc_indices,"hscore"],family="binomial")
decode_gbgc_regression_recomb_hscore=glm(decode_table[gbgc_indices,"fixed|poly"]~decode_table[gbgc_indices,"recombination rate"]*decode_table[gbgc_indices,"hscore"],family="binomial")
decode_gbgc_regression_all=glm(decode_table[gbgc_indices,"fixed|poly"]~decode_table[gbgc_indices,"recombination rate"]+factor(decode_table[gbgc_indices,"genic|nongenic"])+decode_table[gbgc_indices,"hscore"]+decode_table[gbgc_indices,"recombination rate"]*decode_table[gbgc_indices,"hscore"],family="binomial")

summary(decode_gbgc_regression_recomb)
summary(decode_gbgc_regression_genic)
summary(decode_gbgc_regression_hscore)
summary(decode_gbgc_regression_recomb_hscore)
summary(decode_gbgc_regression_all)

save(decode_gbgc_regression_recomb,decode_gbgc_regression_genic,decode_gbgc_regression_hscore,decode_gbgc_regression_recomb_hscore,decode_gbgc_regression_all,file="/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/rdata/decode_recomb_regression_gbgc.Rdata")

##Logistic regressions with HAPMAP recombination rates
##hapmap_table -> HAPMAP recombination rates
hapmap_table=read.table("/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/processed_data/hapmap_recomb_fixed-poly_genic-nongenic_score.bed")
colnames(hapmap_table)=c("chr","start pos.","stop pos.","mut.pattern","recombination rate","fixed|poly","genic|nongenic","hscore")
hapmap_table=hapmap_table[,4:8]
head(hapmap_table)

##all mutation types
##factors separate and with interaction
hapmap_regression_recomb=glm(hapmap_table[,"fixed|poly"]~hapmap_table[,"recombination rate"],family="binomial")
hapmap_regression_genic=glm(hapmap_table[,"fixed|poly"]~factor(hapmap_table[,"genic|nongenic"]),family="binomial")
hapmap_regression_hscore=glm(hapmap_table[,"fixed|poly"]~hapmap_table[,"hscore"],family="binomial")
hapmap_regression_recomb_hscore=glm(hapmap_table[,"fixed|poly"]~hapmap_table[,"recombination rate"]*hapmap_table[,"hscore"],family="binomial")
hapmap_regression_all=glm(hapmap_table[,"fixed|poly"]~hapmap_table[,"recombination rate"]+factor(hapmap_table[,"genic|nongenic"])+hapmap_table[,"hscore"]+hapmap_table[,"recombination rate"]*hapmap_table[,"hscore"],family="binomial")

summary(hapmap_regression_recomb)
summary(hapmap_regression_genic)
summary(hapmap_regression_hscore)
summary(hapmap_regression_recomb_hscore)
summary(hapmap_regression_all)

save(hapmap_regression_recomb,hapmap_regression_genic,hapmap_regression_hscore,hapmap_regression_recomb_hscore,hapmap_regression_all,file="/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/rdata/hapmap_recomb_regression_all.Rdata")

##betwen C<->G
cg_indices_hapmap=grep("[A-Z]{2}[C|G][A-Z]{2}->[A-Z]{2}[C|G][A-Z]{2}",hapmap_table[,1])

hapmap_cg_regression_recomb=glm(hapmap_table[cg_indices_hapmap,"fixed|poly"]~hapmap_table[cg_indices_hapmap,"recombination rate"],family="binomial")
hapmap_cg_regression_genic=glm(hapmap_table[cg_indices_hapmap,"fixed|poly"]~factor(hapmap_table[cg_indices_hapmap,"genic|nongenic"]),family="binomial")
hapmap_cg_regression_hscore=glm(hapmap_table[cg_indices_hapmap,"fixed|poly"]~hapmap_table[cg_indices_hapmap,"hscore"],family="binomial")
hapmap_cg_regression_recomb_hscore=glm(hapmap_table[cg_indices_hapmap,"fixed|poly"]~hapmap_table[cg_indices_hapmap,"recombination rate"]*hapmap_table[cg_indices_hapmap,"hscore"],family="binomial")
hapmap_cg_regression_all=glm(hapmap_table[cg_indices_hapmap,"fixed|poly"]~hapmap_table[cg_indices_hapmap,"recombination rate"]+factor(hapmap_table[cg_indices_hapmap,"genic|nongenic"])+hapmap_table[cg_indices_hapmap,"hscore"]+hapmap_table[cg_indices_hapmap,"recombination rate"]*hapmap_table[cg_indices_hapmap,"hscore"],family="binomial")

summary(hapmap_cg_regression_recomb)
summary(hapmap_cg_regression_genic)
summary(hapmap_cg_regression_hscore)
summary(hapmap_cg_regression_recomb_hscore)
summary(hapmap_cg_regression_all)

save(hapmap_cg_regression_recomb,hapmap_cg_regression_genic,hapmap_cg_regression_hscore,hapmap_cg_regression_recomb_hscore,hapmap_cg_regression_all,file="/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/rdata/hapmap_recomb_regression_cg.Rdata")

##between AT->CG patterns (gBGC control)
gbgc_indices_hapmap=grep("[A-Z]{2}[A|T][A-Z]{2}->[A-Z]{2}[C|G][A-Z]{2}", hapmap_table[,1])

hapmap_gbgc_regression_recomb=glm(hapmap_table[gbgc_indices_hapmap,"fixed|poly"]~hapmap_table[gbgc_indices_hapmap,"recombination rate"],family="binomial")
hapmap_gbgc_regression_genic=glm(hapmap_table[gbgc_indices_hapmap,"fixed|poly"]~factor(hapmap_table[gbgc_indices_hapmap,"genic|nongenic"]),family="binomial")
hapmap_gbgc_regression_hscore=glm(hapmap_table[gbgc_indices_hapmap,"fixed|poly"]~hapmap_table[gbgc_indices_hapmap,"hscore"],family="binomial")
hapmap_gbgc_regression_recomb_hscore=glm(hapmap_table[gbgc_indices_hapmap,"fixed|poly"]~hapmap_table[gbgc_indices_hapmap,"recombination rate"]*hapmap_table[gbgc_indices_hapmap,"hscore"],family="binomial")
hapmap_gbgc_regression_all=glm(hapmap_table[gbgc_indices_hapmap,"fixed|poly"]~hapmap_table[gbgc_indices_hapmap,"recombination rate"]+factor(hapmap_table[gbgc_indices_hapmap,"genic|nongenic"])+hapmap_table[gbgc_indices_hapmap,"hscore"]+hapmap_table[gbgc_indices_hapmap,"recombination rate"]*hapmap_table[gbgc_indices_hapmap,"hscore"],family="binomial")

summary(hapmap_gbgc_regression_recomb)
summary(hapmap_gbgc_regression_genic)
summary(hapmap_gbgc_regression_hscore)
summary(hapmap_gbgc_regression_recomb_hscore)
summary(hapmap_gbgc_regression_all)

save(hapmap_gbgc_regression_recomb,hapmap_gbgc_regression_genic,hapmap_gbgc_regression_hscore,hapmap_gbgc_regression_recomb_hscore,hapmap_gbgc_regression_all,file="/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/rdata/hapmap_recomb_regression_gbgc.Rdata")

