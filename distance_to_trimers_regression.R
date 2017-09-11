##2017-07-26
##regression of distance to closest trimer effect on fixation
##input: patterns_with_closest_distance_all.txt (result of python script "find_distance.py")
##output: 2017_07-26-sorted_regression_distance_to_trimers.txt
##all sorted version

distance_to_trimers_all=read.table("/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/processed_data/patterns_with_closest_distance_all.txt", sep="\t")
colnames(distance_to_trimers_all)=c("chr","start","stop","pattern","rec.rate","fixed|poly","genic|nongenic","hscore","distance")
head(distance_to_trimers_all)

#remove between -3 and 5 (to remove trimers in the pattern itself)
distance_to_trimers_all_2=distance_to_trimers_all[-which((distance_to_trimers_all[,9]>-3) & ((distance_to_trimers_all[,9]<5))),]
#convert negative ones to positive 
distance_to_trimers_all_2[which(distance_to_trimers_all_2[,9]<0),9]=distance_to_trimers_all_2[which(distance_to_trimers_all_2[,9]<0),9]*-1
head(distance_to_trimers_all_2)

#regression
distance_regression=glm(distance_to_trimers_all_2$`fixed|poly`~log(distance_to_trimers_all_2$distance),family=binomial())
summary(distance_regression)

#distribution of distances
hist(distance_to_trimers_all_2[,9])
hist(log(distance_to_trimers_all_2$distance))
hist(distance_to_trimers_all_2[,9])
sum(distance_to_trimers_all_2[,9]<0)
max(distance_to_trimers_all_2[,9])

#correlation
cor.test(distance_to_trimers_all_2$`fixed|poly`,distance_to_trimers_all_2$distance,method="spearman")
