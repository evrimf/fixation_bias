##2017-07-09
##correlation calculation with C<->G patterns and OR>=1
##quint_dfs_human.Rdata includes quint_df_ancestral and quint_df_derived tables


load("/Users/evrimfer/Desktop/evrims/compevo/projeler/Mutation_Fixation_Bias/rdata/quint_dfs_human.Rdata")
head(quint_df_derived)
mut = paste(
  substr(rownames(quint_df_derived), 3, 3), 
  substr(rownames(quint_df_derived), 10, 10) )
table(mut)
x = data.frame( quint_table_derived[ mut == "C G", ] )
par(mfrow=c(1,1))
plot( x$qc1.qc2, x$OR, log="xy", col = grey(1 - (x$score + 3)/7), pch=19)
cor.test(x$qc1.qc2, x$OR, method = "s")
xx = x[x$OR >=1, ]
dim(xx)
cor.test(xx$qc1.qc2, xx$OR, method = "s")

