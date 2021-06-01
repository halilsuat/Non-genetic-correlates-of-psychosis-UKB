library(corrplot)
library(polycor)
library(stringr)

windowsFonts(Times=windowsFont("Times New Roman"))

results <- data.table::fread("results.csv",data.table = F,stringsAsFactors = F)
sigvariables<-results$Predictor[which(results$sig_in_both_datasets==1)]
df_ready_dataset <- data.table::fread("final_Prepared_Dataset.csv",data.table = F,stringsAsFactors = F)
df_outcome_covariates <- data.table::fread("PE_and_covariates.csv",data.table = F,stringsAsFactors = F)
df<-merge(df_ready_dataset,df_outcome_covariates,by="eid")
df<-subset(df,select=-c(eid))
dfchar<- df %>% dplyr::select(where(is.character))
df[names(dfchar)] <- df[names(dfchar)] %>% dplyr::mutate_all(as.factor)
df$sex_f31_0_0<-as.factor(df$sex_f31_0_0)
df$PE<-as.factor(df$PE)
remove(dfchar)
final_multivariate_variables<-c(outcome_and_covariate,sigvariables)
outcome_and_covariate<-c("PE","sex_f31_0_0","age_at_recruitment_f21022_0_0")
final_multivariate_variables<-c(outcome_and_covariate,sigvariables)
df<-df[,grep(paste(final_multivariate_variables, collapse="|"),names(df))]
df<- df[complete.cases(df), ]

df<-subset(df,select=-c(PE))

final_variables_corr_model<-hetcor(df,std.err=F,use="pairwise.complete.obs",pd=F)
final_variables_corr_matric<-final_variables_corr_model$correlations

corrplot()
corr

plot.new()

heatmap.2(final_variables_corr_matric)


final_results <- data.table::fread("final_multivariateresults.csv",data.table = F,stringsAsFactors = F)

final_sig_variables<-final_results$Predictor[which(final_results$Final_significance==1)]

final_sig_variables<-lapply(final_sig_variables,function(x)str_replace_all(x,fixed("Yes"),""))
final_sig_variables<-unlist(final_sig_variables)
final_sig_variables<-final_sig_variables[-1]

df<-subset(df,select=c(final_sig_variables))

final_sig_variables_corr_model<-hetcor(df,std.err=F,use="pairwise.complete.obs",pd=F)
final_sig_variables_corr_matric<-final_sig_variables_corr_model$correlations

heatmapColors<-function(numColors=36){
  c1 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=4/6,end=4.0001/6);
  c2 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=1/6,end=1.0001/6);
  c3 <- c(c1,rev(c2));
  return(c3)
}


pdf("replicated_heat_map.pdf",height=10,width=10)
replicated_heat_map<-heatmap.2(final_sig_variables_corr_matric,trace="none",col=heatmapColors(36),main="replicated",margins=c(36,36))
replicated_heat_map<-replicated_heat_map + theme(text=element_text(family="Times",size=12))
dev.off()

pdf("sig_in_discovery_and_replication_heat_map.pdf",height=10,width=10)
replicated_heat_map<-heatmap.2(final_variables_corr_matric,trace="none",col=heatmapColors(36),main="replicated",margins=c(36,36))
replicated_heat_map<-replicated_heat_map + theme(text=element_text(family="Times",size=12))
dev.off()


corrplot(final_sig_variables_corr_matric, method="circle",type="upper")
corr
