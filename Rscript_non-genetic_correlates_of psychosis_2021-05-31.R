###Author: HSSaraç, 
###Date: 31.5.2021
###Script for detecting the non-genetic correlates of Psychosis. See UK Biobank open-public data showcase for the non-genetic variables

#df is prefered as a generic name for dataframe

#Index of steps:
  
#24 #Step 1: Use data showcase to remove the unuseful Itemtypes, Valuetypes, Strata, and sexed variables + Keep only the first array
#140 #Step 2: Create the outcome variable 
#161 #Step 3: Remove the participants who has a missing value for the outcome variable. 
#171 #Step 4: Uniting of instances
#322 #Step 5: Missing rate cut
#342 #Step 6: Recoding values for each varaible category
#507 #Step 7: Randomly split the dataset into train and test datasets
#535 #Step 8: Logistic Regression
#704 #Step 9: Plotting
  
#Get the packages required
install.packages(c("data.table","bit64","stringr","polycor","DescTools","purrr","dplyr","tidyr","caret","blorr","rcompanion","writexl","ggplot2","ggrepel"))
x<-c("data.table","bit64","stringr","polycor","DescTools","purrr","dplyr","tidyr","caret","blorr","rcompanion","writexl","ggplot2","ggrepel")
lapply(x,require, character.only=TRUE)

#Step 1: Use data showcase to remove the unuseful Itemtypes, Valuetypes, Strata, and sexed variables
df <- data.table::fread("Data_Dictionary_Showcase.csv",data.table = F,stringsAsFactors = F)

#Make the first row: the column names. 
names(df) <- df %>% slice(1) %>% unlist()
df <- df %>% slice(-1)

#Change the column names to how they are named when the dataset is uploaded to the R environment.
df<-data.frame(lapply(df,function(x)tolower(x)))
df<-data.frame(lapply(df,function(x)gsub("[ ]","_",x)))
rownames(df)<-str_c(df$Field,df$FieldID,sep="_f")

patterns=c("-_",",","-","(",")","[","]","/",":","'")
for (pattern in patterns){rownames(df)<-data.frame(lapply(rownames(df),function(x)str_replace_all(x,fixed(pattern),"")))}
rownames(df) <- data.frame(lapply(rownames(df),function(x)str_replace_all(x,fixed("__"),"_")))

#Keep a dataframe with all the variables for later to compare with the removed dataset. Aslýnda sadece rownames tutabilirim.
all_columns<-rownames(df)

Strata_to_remove<- c("auxiliary","supporting") 
df<- df[-which(df$Strata %in% c(Strata_to_remove)),]
Itemtype_to_remove<-c("bulk","samples","records")
df<- df[-which(df$ItemType %in% c(Itemtype_to_remove)),]
Valuetype_to_remove<-c("compound","date","time","binary_object","records")
df<- df[-which(df$ValueType %in% c(Valuetype_to_remove)),]
Sexed_to_remove<-c("female","male")
df<- df[-which(df$Sexed %in% c(Sexed_to_remove)),]

#filter follow-up questions
df_followup<- df %>% 
  filter(str_detect(Notes, "except_those|asked_only|were_asked_if_the_participant_indicated|only_asked|as_defined_by_their_answer|was_only_asked_to_participants_who|question_was_asked_when"))
#remove these questions by filtering out the rows 
df<-df[!(rownames(df) %in% rownames(df_followup)), ]

#the list of variables to be removed.
remove<- all_columns [!(all_columns%in% rownames(df))] 

#Turn array variable into numeric class.
df$Array<-as.numeric(df$Array)
#Retrieve the variables with multiple arrays
multiple_array_items_only <- which(df$Array>1)
dfarray <- df[multiple_array_items_only, ]
Items_with_multiple_arrays<-rownames(dfarray)

#Retrieve the numeric variables
numeric_items_only <- which(df$ValueType=="continuous"|df$ValueType=="integer")
df <- df[numeric_items_only, ]
Items_with_numeric<-rownames(df)

remove <- as.data.frame(remove)

Items_with_multiple_arrays <- as.data.frame(Items_with_multiple_arrays)
Items_with_numeric <- as.data.frame(Items_with_numeric)

remove(df,numeric_items_only,dfarray)

#Opening the ukb_data.csv file on R as a dataframe using fread

df <- data.table::fread("tinyukb.csv",data.table = F,stringsAsFactors = F)

namesdf<-names(df)
write.csv(namesdf,file="allnames_df.csv")
remove(namesdf)

#remove the auxiliary, bulk, compound, follow-up, time, date and raw data variables.
for(i in 1:nrow(remove)){
  removea<-remove[i,]
  df<-df[,!grepl(paste(removea,collapse="|"),names(df))]
}
#Remove other auxiliary and genetic variables
removed_auxiliary_genetic<- c("affymetrix_quality_control_metric_clustercr_f22025","affymetrix_quality_control_metric_dqc_f22026","genotype_measurement_batch_f22000","genetic_kinship_to_other_participants_f22021","genetic_sex_f22001","genotype_measurement_plate_f22007","genotype_measurement_well_f22008","index_for_card_a_in_round_f401","index_for_card_b_in_round_f402","month_of_attending_assessment_centre_f55","ordering_of_blows_f3065","spirometry_method_f23","uk_biobank_assessment_centre_f54","was_blood_sampling_attempted_f35","weight_method_f21","age_when_attended_assessment_centre_f21003","fasting_time_f74","genetic_principal_components_f22009","heterozygosity_f22003","heterozygosity_pca_corrected_f22004","home_location_at_assessment_east_coordinate_rounded_f20074","home_location_at_assessment_north_coordinate_rounded_f20075","home_location_east_coordinate_rounded_f22702","home_location_north_coordinate_rounded_f22704","missingness_f22005","number_of_blood_samples_taken_f68","number_of_measurements_made_f3137","number_of_times_snapbutton_pressed_f403","place_of_birth_in_uk_east_coordinate_f130","place_of_birth_in_uk_north_coordinate_f129","time_to_complete_round_f400","year_of_birth_f34","number_of_correct_matches_in_round_f398","birth_weight_known_f120","V1")
df=df[,!grepl(paste(removed_auxiliary_genetic, collapse = "|"), names(df))]
remove(removed_auxiliary_genetic)
remove(df_followup,all_columns,multiple_array_items_only,Sexed_to_remove,Strata_to_remove,Valuetype_to_remove,patterns,pattern,i,Itemtype_to_remove)
remove(remove,removea)

#Only keep variables that are the first array of a multi array instance. The second numeric value in the list below represents the array number, 0 is the first array. The list ends with "_31_0", because library(sjmisc) find_var(df,"_32_0") gives an empty argument suggesting there is no such variable. 
#make Items with arrays the rowname , #job_soc_coding is added because it is named differently as "job code historical" in the showcase.
Items_with_multiple_arrays<-Items_with_multiple_arrays$Items_with_multiple_arrays
Items_with_multiple_arrays<-c(Items_with_multiple_arrays,"job_soc_coding_f22617")

#Dataset of multi array variables
dfmultiple_arrays<-df[,grep(paste(c(Items_with_multiple_arrays,"eid"), collapse="|"),names(df))]
firstarrays<-c("_0_0","_1_0","_2_0","_3_0","_4_0","_5_0","_6_0","_7_0","_8_0","_9_0","_10_0","_11_0","_12_0","_13_0","_14_0","_15_0","_16_0","_17_0","_18_0","_19_0","_20_0","_21_0","_22_0","_23_0","_24_0","_25_0","_26_0","_27_0","_28_0","_29_0","_30_0","_31_0","eid")
#only keep the first array
df_multiple_arrays_normal<-dfmultiple_arrays[,grep(paste(firstarrays, collapse="|"),names(dfmultiple_arrays))]
remove(firstarrays)
Items_with_multiple_arrays_normal<-names(df_multiple_arrays_normal)
#Delete _0_0 instance array codes so that we only have the variable name and fcode.
patterns=c("_0_0","_1_0","_2_0","_3_0")
for (pattern in patterns){Items_with_multiple_arrays_normal<-(lapply(Items_with_multiple_arrays_normal,function(x)str_replace_all(x,fixed(pattern),"")))}

Items_with_deviant_multiple_arrays <- Items_with_multiple_arrays[!(Items_with_multiple_arrays%in% Items_with_multiple_arrays_normal)] 
df_multiple_arrays_deviant<-dfmultiple_arrays[,grep(paste(c(Items_with_deviant_multiple_arrays,"eid"), collapse="|"),names(dfmultiple_arrays))]

firstarrays<-c("_0_1","_1_1","_2_1","_3_1","eid")
#only keep the first array
df_multiple_arrays_deviant<-df_multiple_arrays_deviant[,grep(paste(firstarrays, collapse="|"),names(df_multiple_arrays_deviant))]
todelete<-c("_0_10","_0_11","_0_12","_0_13","_0_14","_0_15","_0_16","_0_17","_0_18","_1_10","_1_11","_1_12","_1_13","_1_14","_1_15","_1_16","_1_17","_1_18","_2_10","_2_11","_2_12","_2_13","_2_14","_2_15","_2_16","_2_17","_2_18","_3_10","_3_11","_3_12","_3_13","_3_14","_3_15","_3_16","_3_17","_3_18")
df_multiple_arrays_deviant<-df_multiple_arrays_deviant[,!grepl(paste(todelete, collapse="|"),names(df_multiple_arrays_deviant))]

#Create a dataset with one array variables
df_single_array_variables<-df[,!grepl(paste(Items_with_multiple_arrays, collapse="|"),names(df))]
remove(df)

dfmultiple_arrays<-merge(df_multiple_arrays_normal,df_multiple_arrays_deviant,by="eid")
remove(df_multiple_arrays_normal,df_multiple_arrays_deviant,todelete,firstarrays)

df<-merge(df_single_array_variables,dfmultiple_arrays,by="eid")
remove(df_single_array_variables,dfmultiple_arrays,Items_with_multiple_arrays,Items_with_deviant_multiple_arrays,Items_with_multiple_arrays_normal)

#tarpd variables and source variables are renamed to avoide confusion with "f1","f2" when uniting instances _f denotes to the start of the field ID.
names(df) <- gsub("tarpd_f1_antigen_for_chlamydia_trachomatis_f23035_0_0", "tarpdf1_antigen_for_chlamydia_trachomatis_f23035_0_0",names(df))
names(df) <- gsub("tarpd_f2_antigen_for_chlamydia_trachomatis_f23036_0_0", "tarpdf2_antigen_for_chlamydia_trachomatis_f23036_0_0",names(df))
names(df) <- gsub("source_of_report_of_f", "source_of_report_off",names(df))

##Step 2: Create the outcome variable uniting of psychotic experiences(PE) variables
#Values of "Prefer not to answer" and "Do not know" are coded as NA across the PE variables. 
df$ever_heard_an_unreal_voice_f20463_0_0[df$ever_heard_an_unreal_voice_f20463_0_0 == "Do not know" | df$ever_heard_an_unreal_voice_f20463_0_0 == "Prefer not to answer"]<-NA
df$ever_seen_an_unreal_vision_f20471_0_0[df$ever_seen_an_unreal_vision_f20471_0_0 == "Do not know" | df$ever_seen_an_unreal_vision_f20471_0_0 == "Prefer not to answer"]<-NA
df$ever_believed_in_an_unreal_conspiracy_against_self_f20468_0_0[df$ever_believed_in_an_unreal_conspiracy_against_self_f20468_0_0 == "Do not know" | df$ever_believed_in_an_unreal_conspiracy_against_self_f20468_0_0 == "Prefer not to answer"]<-NA
df$ever_believed_in_unreal_communications_or_signs_f20474_0_0[df$ever_believed_in_unreal_communications_or_signs_f20474_0_0 == "Do not know" | df$ever_believed_in_unreal_communications_or_signs_f20474_0_0 == "Prefer not to answer"]<-NA

#PE variables are united into a single column:P. NA is infectious in stringr, any column that has even a single NA column turns into NA.
df$P <- stringr::str_c(df$ever_heard_an_unreal_voice_f20463_0_0,df$ever_seen_an_unreal_vision_f20471_0_0,df$ever_believed_in_an_unreal_conspiracy_against_self_f20468_0_0,df$ever_believed_in_unreal_communications_or_signs_f20474_0_0)

##Creating the outcome column: "PE": Coding as "0" the noPE group who responded "No" to all of the PE questions.
df$PE<- ifelse(df$P == "NoNoNoNo", "No", "Yes")

##If any of the 4 PE questions is answered "Yes" we coded PE as "1".So, we have also included rows which included both "Yes" and NA(which was not accounted for in stringr for which NA is infectious).
df$PE[df$ever_heard_an_unreal_voice_f20463_0_0=="Yes"|df$ever_seen_an_unreal_vision_f20471_0_0=="Yes"|df$ever_believed_in_an_unreal_conspiracy_against_self_f20468_0_0=="Yes"|df$ever_believed_in_unreal_communications_or_signs_f20474_0_0=="Yes"]<-"Yes"

#After this analysis, we end up with 7803 yesPE and 147444 noPE participants, the rest is NA.

##Delete the strings used to create the outcome variable from the dataset.
df = subset(df, select = -c(ever_heard_an_unreal_voice_f20463_0_0,ever_seen_an_unreal_vision_f20471_0_0,ever_believed_in_an_unreal_conspiracy_against_self_f20468_0_0,ever_believed_in_unreal_communications_or_signs_f20474_0_0,P) )

###Step 3: Delete the rows that has a missing value for the outcome variable. 
df<- df[!is.na(df$PE),]

#Making a dataset of only the participant id, PE, and covariates
df$PE<- as.factor(df$PE)
PE = subset(df, select = c("PE","eid","age_at_recruitment_f21022_0_0","sex_f31_0_0") )
write.csv(PE,file="PE_and_covariates.csv",row.names=FALSE)
remove(PE)
df<-subset(df,select=-c(PE,age_at_recruitment_f21022_0_0,sex_f31_0_0))

#Step 4: Uniting of instances 
#seperate the dataset into character variables and numerical variables
dfchar<-df
for(i in 1:nrow(Items_with_numeric)){
  Items_with_numerica<-Items_with_numeric[i,]
  dfchar<-dfchar[,!grepl(paste(Items_with_numerica,collapse="|"),names(dfchar))]
}

dfchar <- dfchar %>% dplyr::mutate_all(as.character)
namesdfnum<- names(df)[!(names(df)%in% names(dfchar))] 
eid<-"eid"
namesdfnum<-c(eid,namesdfnum)
dfnum<-subset(df,select=c(namesdfnum))

remove(df,Items_with_numeric,Items_with_numeric1,Items_with_numeric2,Items_with_numeric3,Items_with_numeric4,Items_with_numeric5,Items_with_numeric6,Items_with_numeric7,Items_with_numeric8,Items_with_numeric9,Items_with_numeric10,Items_with_numeric11,Items_with_numeric12,Items_with_numeric13,Items_with_numeric14,Items_with_numeric15,Items_with_numeric16,Items_with_numeric17,eid,namesdfnum)

# for loop dflist could also be used for this part. 

# Uniting instances for the character variables
#Leave no space in variable values, so that they are recognized as one value.
dfchar<-data.frame(lapply(dfchar,function(x)gsub("[ ]","_",x)))
#Create a space between variable name and f-code
names(dfchar)<-gsub("\\_(_*?)(?=f\\d)"," ",names(dfchar),perl=T)

#Change the variable names in a way that the different instances of the same variable are grouped together when uniting columns.
#Split the name of the column by space which seperates name of the variable(V1) from the fcode + instance and array number(V2) of the name 
Names <- as.data.frame(str_split_fixed(colnames(dfchar[-1]), " ", 2))
#Turn the "_" into space
Names$V2 <- gsub("_", " ",Names$V2)
#Split the f code and the instance array numbers
Names_field_id <- as.data.frame(str_split_fixed(Names$V2, " ", 2))
#keep the f code 
Names$V2 <- word(Names$V2)
#Unite variable name with the f code, this(Var. name+f-code) is the pattern we wil be using to unite instances.
Names$V3 <- stringr::str_c(Names$V1,Names$V2,sep="_")
#Also add the instance and array numbers seperated by space
Names$V4 <- stringr::str_c(Names$V3,Names_field_id$V2, sep=" ")
colnames(dfchar) <- c("eid",Names$V4)

#Now we can create cols where the variables with different instances are clustered.
cols <- split(names(dfchar)[-1],word(names(dfchar)[-1], 1))
#The instances of a variable are binded together, after each value, there is a space so that we can later seperate instances. 
dfchar<- cbind(dfchar[1],sapply(names(cols), function(col) {
  do.call(paste, c(dfchar[cols[[col]]], sep = " ",na.rm=TRUE))
}))

#NA and TRUE are removed.
remove(cols)
patterns=c(" NA","NA "," TRUE","TRUE ","TRUE")
for (pattern in patterns){dfchar<-data.frame(lapply(dfchar,function(x)str_replace_all(x,fixed(pattern),"")))}

eid<-c("eid")
dfcharvars<-dfchar[,!grepl(paste(eid, collapse="|"),names(dfchar))]
#Only the first value is kept.
dfchar[names(dfcharvars)]<-data.frame(lapply(dfchar[names(dfcharvars)],function(x)word(x,1)))

#When there is no value at all, it is recoded as space, this space is deleted.
dfchar<-data.frame(lapply(dfchar,function(x)gsub("[ ]","",x)))
#Variables with NA in all instances are coded as NA
dfchar[dfchar==""]<-NA

remove(dfcharvars,eid,Names,Names_field_id)

#Leave no space in variable values, so that they are recognized as one value.
dfnum<-data.frame(lapply(dfnum,function(x)gsub("[ ]","_",x)))
#Create a space between variable name and f-code
names(dfnum)<-gsub("\\_(_*?)(?=f\\d)"," ",names(dfnum),perl=T)

Names <- as.data.frame(str_split_fixed(colnames(dfnum[-1]), " ", 2))
Names$V2 <- gsub("_", " ",Names$V2)
Names_field_id <- as.data.frame(str_split_fixed(Names$V2, " ", 2))
Names$V2 <- word(Names$V2)
Names$V3 <- stringr::str_c(Names$V1,Names$V2,sep="_")
Names$V4 <- stringr::str_c(Names$V3,Names_field_id$V2, sep=" ")
colnames(dfnum) <- c("eid",Names$V4)

cols <- split(names(dfnum)[-1],word(names(dfnum)[-1], 1))

dfnum<- cbind(dfnum[1],sapply(names(cols), function(col) {
  do.call(paste, c(dfnum[cols[[col]]], sep = " ",na.rm=TRUE))
}))

remove(cols)
patterns=c(" NA","NA "," TRUE","TRUE ","TRUE")
for (pattern in patterns){dfnum<-data.frame(lapply(dfnum,function(x)str_replace_all(x,fixed(pattern),"")))}

eid<-c("eid")
dfnumvars<-dfnum[,!grepl(paste(eid, collapse="|"),names(dfnum))]

dfnum[names(dfnumvars)]<-data.frame(lapply(dfnum[names(dfnumvars)],function(x)word(x,1)))

dfnum<-data.frame(lapply(dfnum,function(x)gsub("[ ]","_",x)))

dfnum[dfnum==""]<-NA

remove(dfnumvars,eid,Names,Names_field_id)

dfnum <- dfnum %>% dplyr::mutate_all(as.numeric)
dfchar <- dfchar %>% dplyr::mutate_all(as.character)

#merge the dataset back.
df<-merge(dfchar,dfnum,by="eid")

#Recoding "_" as space back to the original values of the dataset.
df<-data.frame(lapply(df,function(x)gsub("[_]"," ",x)))

df[names(dfnum)] <- df[names(dfnum)] %>% dplyr::mutate_all(as.numeric)
df[names(dfchar)] <- df[names(dfchar)] %>% dplyr::mutate_all(as.character)

remove(dfchar,dfnum)

#Variables with no variance are removed. 
dffactorized<-df %>% dplyr::mutate_all(as.factor)
nlevels<-lapply(dffactorized,nlevels)

#Remove all variables with single level
single_or_no_level_variables<-nlevels[nlevels==1]
single_or_no_level_variables<-as.data.frame(single_or_no_level_variables)
df<-df[,!(names(df) %in% names(single_or_no_level_variables)) ]

write.csv(single_or_no_level_variables,file="singlelevels.csv",row.names = FALSE)

remove(dffactorized,single_or_no_level_variables,nlevels)

df[df == "Do not know" | df == "Prefer not to answer"|df == "Unsure"] <- NA
df[df=="Do not drive on the motorway"|df=="Do not drink hot drinks"|df=="Post code not linkable"] <- NA

#NA coding of numerical variables
dfnum<-df %>% dplyr::select(where(is.numeric))

df_numerical_factorized<-df[names(dfnum)] %>% dplyr::mutate_all(as.factor)
values<-lapply(df_numerical_factorized,levels)
numberofneg<-lapply (values,function(x)length(which(x<0)))
highnegs<-numberofneg[numberofneg>3]
numeric_with_negative<-c(names(highnegs))
lownegs<-numberofneg[numberofneg<=3]
numeric_with_positive<-c(names(lownegs))

remove(df_numerical_factorized,dfnum)
#Recoding of values #Recode placeholder values as NA based on the coding showcase. Maybe I can write their data-coding number by them.
df[numeric_with_positive][df[numeric_with_positive]=="-1"|df[numeric_with_positive]=="-2"|df[numeric_with_positive]=="-3"]<-NA

remove(numeric_with_positive,highnegs,lownegs,numberofneg,values,numeric_with_negative)

#missing rate calculation
missing_rate1_list_num <- sapply(df,function(x) sum(is.na(x)/nrow(df)))
dfmissing<-as.data.frame(missing_rate1_list_num)
dfmissing$missing_rate1_list_num<-dfmissing$missing_rate1_list_num*100

write.csv(dfmissing,file="missingratetable.csv")

#Step 5: Missing Rate Cut
df<-df %>%  purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=10)
remove(dfmissing,missing_rate1_list_num)
write.csv(df,file="variables_after_missing_cut.csv",row.names=FALSE)

df_predictors <- data.table::fread("Predictor-table.csv",data.table = F,stringsAsFactors = F)
rownames(df_predictors)<-str_c(df_predictors$`Predictor Variable`,df_predictors$`f code`, sep="_f")

#Make the first row: the column names. 
find_rows <- which(df_predictors$`Variable Type` == "Categorical Ordered")
df_ordered_categorical <- df_predictors[find_rows, ]

find_rows <- which(df_predictors$`Variable Type` == "Categorical Unordered")
df_unordered_categorical <- df_predictors[find_rows, ]

find_rows <- which(df_predictors$`Variable Type` == "Numerical")
df_numerical <- df_predictors[find_rows, ]

remove(find_rows,df_predictors)

#Step 6: Recoding values for each varaible category
df$eid<-as.character(df$eid)

#a) Un-ordered categorical variables recoding
# Value coding: It will be recoded as "No" for the most populous group and "Yes" for the rest. 

#Detecting the most populated value in a variable

populous_groups<-as.data.frame(lapply(df[rownames(df_unordered_categorical)],function(x)Mode(x,na.rm=TRUE)))
write.csv(populous_groups,file="populous_Groups_for_non_ordered_variables.csv")
for(variable in colnames(populous_groups)) 
  df[,variable]= ifelse(df[,variable] == populous_groups[1,variable],"No","Yes") 
remove(populous_groups)

#b) coding of the Ordered categorical variables
#Ordered categorical variables are recoded. The lowest value of the variable:0, the rest:1

df_predictors <- data.table::fread("Predictor-table.csv",data.table = F,stringsAsFactors = F)
rownames(df_predictors)<-str_c(df_predictors$`Predictor Variable`,df_predictors$`f code`, sep="_f")
find_rows <- which(df_predictors$`Variable Type` == "Categorical Ordered")
df_predictors <- df_predictors[find_rows, ]
predictors_table<-subset(df_predictors,select="Recoded Values")
remove(df_predictors)
#prepare the predictor table
patterns=c(" |","| ")
for (pattern in patterns){predictors_table<-data.frame(lapply(predictors_table,function(x)str_replace_all(x,fixed(pattern),"|")))}

patterns=c(" =","= ")
for (pattern in patterns){predictors_table<-data.frame(lapply(predictors_table,function(x)str_replace_all(x,fixed(pattern),"=")))}
#replace | with __ because | causes problem in the strSplit step. 
predictors_table<-as.data.frame(lapply(predictors_table,function(x)gsub("[|]","__",x)))
values<- StrSplit(predictors_table$Recoded.Values,split="__")
values<- as.data.frame(unlist(values))
#Create old and new values columns by seperating the recoded values by "="
values_table<-as.data.frame(str_split_fixed(values$`unlist(values)`,"=",2))
remove(values)
names(values_table)<-c("old","new")
find_rows_1 <- which(values_table$new == "1")
values_table_1 <- values_table[find_rows_1, ]
find_rows_0 <- which(values_table$new == "0")
values_table_0 <- values_table[find_rows_0, ]
#Take only one value in the case of duplication
values_table_0<-values_table_0%>%group_by(values_table_0$old)%>%distinct
values_table_1<-values_table_1%>%group_by(values_table_1$old)%>%distinct

old_value<-c(values_table_0$old,values_table_1$old)
new_value<-c(values_table_0$new,values_table_1$new)
#For every value change the old one with the new one.
for(i in 1:111){
  df[rownames(df_ordered_categorical)][df[rownames(df_ordered_categorical)]==old_value[i]]<-new_value[i]
}
remove(predictors_table,values_table_0,values_table_1,values_table,new_value,old_value,i,find_rows_0,find_rows,find_rows_1)

#Turn all recodings to "Yes" and "No" so R recognizes the variables as character. 
df[rownames(df_ordered_categorical)][df[rownames(df_ordered_categorical)]=="0"]<-"No"
df[rownames(df_ordered_categorical)][df[rownames(df_ordered_categorical)]=="1"]<-"Yes"

#c)Recoding for numerical
#Seperate the numerical class variables into 2: ones that have (-) values and one that only has (+) values.  This is to avoid unintentionally recoding meaningful (-) values 
df_numerical_factorized<-df[rownames(df_numerical)] %>% dplyr::mutate_all(as.factor)
values<-lapply(df_numerical_factorized,levels)
numberofneg<-lapply (values,function(x)length(which(x<0)))
highnegs<-numberofneg[numberofneg>1]
numeric_with_negative<-c(names(highnegs))
lownegs<-numberofneg[numberofneg<=1]
numeric_with_positive<-c(names(lownegs))

remove(df_numerical_factorized)
#Recoding of values
df[numeric_with_positive]<-df[numeric_with_positive] %>% dplyr::mutate_all(as.character)

# All the codings of -10 are similar they refer to values between 0 and 1(less than one) in integer columns. Most of these variables are numerical though!
df[numeric_with_positive] [df[numeric_with_positive]== "-10" ] <-"0"

df[numeric_with_positive]<-df[numeric_with_positive] %>% dplyr::mutate_all(as.numeric)

#numeric variables with less than 10 levels will be recoded 0-1.
df_numerical_factorized<-df[rownames(df_numerical)] %>% dplyr::mutate_all(as.factor)
nlevels<-lapply(df_numerical_factorized,nlevels)

below10level_variables<-as.data.frame(nlevels[nlevels<10])
minimums_of_less_than_10_values<-lapply(df[colnames(below10level_variables)],min,na.rm=TRUE)
df[rownames(df_numerical)]<-df[rownames(df_numerical)] %>% dplyr::mutate_all(as.numeric)

df[colnames(below10level_variables)] = ifelse(df[colnames(below10level_variables)] =="0", "No","Yes")
df[colnames(below10level_variables)] <- df[colnames(below10level_variables)] %>% dplyr::mutate_all(as.character)

write.csv(minimums_of_less_than_10_values,file="below10level_minimum_values.csv")
remove(minimums_of_less_than_10_values)

write.csv(below10level_variables,file="below10vars.csv")
remove(nlevels,below10level_variables)

#z-score transformation of numeric variables
dfnum<-df %>% dplyr::select(where(is.numeric))
df[names(dfnum)] <-scale(df[names(dfnum)],center=TRUE,scale=TRUE)

remove(dfnum,df_numerical,df_numerical_factorized,df_ordered_categorical,df_unordered_categorical,highnegs,lownegs,numberofneg,values,numeric_with_negative,numeric_with_positive)

write.csv(df,file="PreparedDataset.csv",row.names=FALSE)

#Correlation analysis for detecting highly collinear variables
df_with_eid<-df

#We remove "eid" from the dataset for the collinearity analysis
df<-subset(df,select=-c(eid))

correlation_model<-hetcor(df,std.err=F,use="pairwise.complete.obs",pd=F)
correlation_matric<-correlation_model$correlations
df_correlation<-as.data.frame(correlation_matric)

write.csv(df_correlation,file="Correlation-matrice.csv")

#Find correlations higher than or equal to .90
correlation_high<-ifelse(df_correlation>=.90,1,0)

df_correlation_high<-as.data.frame(correlation_high)

high_cors<-sapply(df_correlation_high,function(x)length(which(x==1)))
high_cors_variables<-high_cors[high_cors>1]

#Find correlations lower than or equal to -.90
correlation_low<-ifelse(df_correlation<=-.90,1,0)
df_correlation_low<-as.data.frame(correlation_low)
low_cors<-sapply(df_correlation_low,function(x)length(which(x==1)))
high_neg_cors_variables<-low_cors[low_cors>0]

#Make a list of variables that are negatively and positively highly(<=-.90,>=.90) correlated
correlated_variables<-c(names(high_cors_variables),names(high_neg_cors_variables))
df_correlated_variables<-as.data.frame(correlated_variables)
corvars<-df_correlated_variables$correlated_variables

#Make a dataset of these highly correlated variables
df_high_cor=df[,grep(paste(corvars, collapse = "|"), names(df))]

remove(df,correlated_variables)

#Correlation matrix of highly correlated variables is saved. 
correlation_matrix_high_cors<-hetcor(df_high_cor,std.err=F,use="pairwise.complete.obs",pd=F)
highcors<-correlation_matrix_high_cors$correlations
highcors<-as.data.frame(highcors)
remove(df_high_cor,correlation_matrix_high_cors)
write_xlsx(highcors,"Supplementary_table_3_highcorvars.xlsx")
remove(highcors)

#Removal of a variable from each highly collinear pair through findCorrelation command that removes variables based on the idea to keep it collinearity of the whole dataset as low as possible. 
highly_collinear_variables_to_remove<-findCorrelation(correlation_matric,cutoff=0.9,names=TRUE)

df_with_eid=df_with_eid[,!grepl(paste(highly_collinear_variables_to_remove, collapse = "|"), names(df_with_eid))]

write.csv(df_with_eid,file="final_Prepared_Dataset.csv",row.names=FALSE)

remove(correlation_high,correlation_low,correlation_matric,correlation_model,df_correlated_variables,df_correlation,df_correlation_high,df_correlation_low,corvars,high_cors,high_cors_variables,high_neg_cors_variables,highly_collinear_variables_to_remove,low_cors)

df<-df_with_eid

df_outcome_covariates <- data.table::fread("PE_and_covariates.csv",data.table = F,stringsAsFactors = F)

df<-merge(df,df_outcome_covariates,by="eid")

#We remove eid from the analysis dataset.
df<-subset(df,select=-c(eid))

remove(df_outcome_covariates)

#Step 7: Randomly split the dataset into train and test datasets
##Split it proportionately to the outcome measure(PE) 
set.seed(123)
discovery.index <- 
  createDataPartition(df$PE,p = 0.5, list = FALSE)
df_discovery <- df[discovery.index, ]
df_replication <- df[-discovery.index, ]

remove(discovery.index,df)

cols_predictors<-names(df_discovery)
cols_predictors<- cols_predictors[!cols_predictors%in% c("PE","age_at_recruitment_f21022_0_0","sex_f31_0_0")] 
cols_predictors<-as.data.frame(cols_predictors)
rownames(cols_predictors)<-cols_predictors$cols_predictors

#Rearranging the classes of character variables, outcome and covariate
df_discoverychar<- df_discovery %>% dplyr::select(where(is.character))

df_discovery[names(df_discoverychar)] <- df_discovery[names(df_discoverychar)] %>% dplyr::mutate_all(as.factor)
df_replication[names(df_discoverychar)] <- df_replication[names(df_discoverychar)] %>% dplyr::mutate_all(as.factor)

df_discovery$sex_f31_0_0<-as.factor(df_discovery$sex_f31_0_0)
df_replication$sex_f31_0_0<-as.factor(df_replication$sex_f31_0_0)
df_discovery$PE<-as.factor(df_discovery$PE)
df_replication$PE<-as.factor(df_replication$PE)

remove(df_discoverychar)

#Step 8: #Logistic regression
#a) In the training dataset, analyze each predictor's association with PE with the covariates: age and sex.

results_discovery <- data.frame(Predictor=rownames(cols_predictors),Discovery_Pval=NA,Discovery_Coefficient=NA,Discovery_SD=NA,Discovery_Oddsratios=NA,Discovery_R2=NA)

for(i in 1:nrow(cols_predictors)){
  modelA <- summary(glm(PE ~ df_discovery[,i]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df_discovery, family = binomial))
  results_discovery[i,2] <- modelA$coefficients[2,4] ##2.line 4.column
  results_discovery[i,3] <- modelA$coefficients[2,1]
  results_discovery[i,4] <- modelA$coefficients[2,2]
}

#To get the pseudo R2 for each model
R2_glms_discovery<-lapply(1:nrow(cols_predictors),function(x)nagelkerke(glm(PE ~ df_discovery[,x]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df_discovery, family = binomial),null=glm(PE ~ 1, binomial, df_discovery[!is.na(df_discovery[,x]),])))
Nagelkerke_Cragg_Uhler_discovery<- lapply(R2_glms_discovery,function(x) x$Pseudo.R.squared.for.model.vs.null[3,1])
Nagelkerke_Cragg_Uhler_discovery<-unlist(Nagelkerke_Cragg_Uhler_discovery)
results_discovery$Discovery_R2<-Nagelkerke_Cragg_Uhler_discovery

results_discovery$Discovery_Oddsratios<- exp(results_discovery$Discovery_Coefficient)
results_discovery$Discovery_significance<-ifelse(results_discovery$Discovery_Pval<=0.05/nrow(cols_predictors),1,0)
results_discovery$Discovery_ORlci <- exp(results_discovery$Discovery_Coefficient + qnorm(c(0.025)) * results_discovery$Discovery_SD)
results_discovery$Discovery_ORhci <- exp(results_discovery$Discovery_Coefficient + qnorm(c(0.975)) * results_discovery$Discovery_SD)

results_discovery<-as.data.frame(results_discovery)
write.csv(results_discovery,file="results_discovery.csv",row.names = FALSE)

#b)In the Validation dataset
results_replication <- data.frame(Predictor=rownames(cols_predictors),Replication_Pval=NA,Replication_Coefficient=NA,Replication_SD=NA,Replication_Oddsratios=NA,Replication_R2=NA)

for(i in 1:nrow(cols_predictors)){
  modelB <- summary(glm(PE ~ df_replication[,i]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df_replication, family = binomial))
  results_replication[i,2] <- modelB$coefficients[2,4] ##2.line 4.column
  results_replication[i,3] <- modelB$coefficients[2,1]
  results_replication[i,4] <- modelB$coefficients[2,2]
}

#extract the Pseudo R2
R2_glms_replication<-lapply(1:nrow(cols_predictors),function(x)nagelkerke(glm(PE ~ df_replication[,x]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df_replication, family = binomial),null=glm(PE ~ 1, binomial, df_replication[!is.na(df_replication[,x]),])))
Nagelkerke_Cragg_Uhler_replication<- lapply(R2_glms_replication,function(x) x$Pseudo.R.squared.for.model.vs.null[3,1])
Nagelkerke_Cragg_Uhler_replication<-unlist(Nagelkerke_Cragg_Uhler_replication)
results_replication$Replication_R2<-Nagelkerke_Cragg_Uhler_replication

results_replication$Replication_Oddsratios<- exp(results_replication$Replication_Coefficient)
results_replication$Replication_significance<-ifelse(results_replication$Replication_Pval<=0.05/nrow(cols_predictors),1,0)
results_replication$Replication_ORlci <- exp(results_replication$Replication_Coefficient + qnorm(c(0.025)) * results_replication$Replication_SD)
results_replication$Replication_ORhci <- exp(results_replication$Replication_Coefficient + qnorm(c(0.975)) * results_replication$Replication_SD)

write.csv(results_replication,file="results_replication.csv",row.names=FALSE)

remove(df_discovery,df_replication,modelA,modelB,R2_glms_discovery,R2_glms_replication,i,Nagelkerke_Cragg_Uhler_discovery,Nagelkerke_Cragg_Uhler_replication)

#Regression Analysis in the combined dataset
df_ready_dataset <- data.table::fread("final_Prepared_Dataset.csv",data.table = F,stringsAsFactors = F)
df_outcome_covariates <- data.table::fread("PE_and_covariates.csv",data.table = F,stringsAsFactors = F)

df<-merge(df_ready_dataset,df_outcome_covariates,by="eid")

df<-subset(df,select=-c(eid))

dfchar<- df %>% dplyr::select(where(is.character))

df[names(dfchar)] <- df[names(dfchar)] %>% dplyr::mutate_all(as.factor)

df$sex_f31_0_0<-as.factor(df$sex_f31_0_0)
df$PE<-as.factor(df$PE)

remove(dfchar)

results_combined <- data.frame(Predictor=rownames(cols_predictors),Combined_Pval=NA,Combined_Coefficient=NA,Combined_SD=NA,Combined_Oddsratios=NA,Combined_R2=NA)

for(i in 1:nrow(cols_predictors)){
  modelC <- summary(glm(PE ~ df[,i]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df, family = binomial))
  results_combined[i,2] <- modelC$coefficients[2,4] ##2.line 4.column
  results_combined[i,3] <- modelC$coefficients[2,1]
  results_combined[i,4] <- modelC$coefficients[2,2]
}

R2_glms_combined_first_130<-lapply(1:130,function(x)nagelkerke(glm(PE ~ df[,x]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df, family = binomial),null=glm(PE ~ 1, binomial, df[!is.na(df[,x]),])))
Nagelkerke_Cragg_Uhler_combined_first_130 <- lapply(R2_glms_combined_first_130,function(x) x$Pseudo.R.squared.for.model.vs.null[3,1])
Nagelkerke_Cragg_Uhler_combined_first_130 <-unlist(Nagelkerke_Cragg_Uhler_combined_first_130)
results_combined[1:130,]$Combined_R2<-Nagelkerke_Cragg_Uhler_combined_first_130

R2_glms_combined_last_130<-lapply(131:nrow(cols_predictors),function(x)nagelkerke(glm(PE ~ df[,x]+age_at_recruitment_f21022_0_0+sex_f31_0_0, data = df, family = binomial),null=glm(PE ~ 1, binomial, df[!is.na(df[,x]),])))
Nagelkerke_Cragg_Uhler_combined_last_130 <- lapply(R2_glms_combined_last_130,function(x) x$Pseudo.R.squared.for.model.vs.null[3,1])
Nagelkerke_Cragg_Uhler_combined_last_130 <-unlist(Nagelkerke_Cragg_Uhler_combined_last_130)
results_combined[131:nrow(cols_predictors),]$Combined_R2<-Nagelkerke_Cragg_Uhler_combined_last_130

results_combined$Combined_Oddsratios<- exp(results_combined$Combined_Coefficient)
results_combined$Combined_significance<-ifelse(results_combined$Combined_Pval<=0.05/nrow(cols_predictors),1,0)
results_combined$Combined_ORlci <- exp(results_combined$Combined_Coefficient + qnorm(c(0.025)) * results_combined$Combined_SD)
results_combined$Combined_ORhci <- exp(results_combined$Combined_Coefficient + qnorm(c(0.975)) * results_combined$Combined_SD)

write.csv(results_combined,file="results_combined.csv",row.names=FALSE)

remove(modelC,i,Nagelkerke_Cragg_Uhler_combined_first_130,Nagelkerke_Cragg_Uhler_combined_last_130,R2_glms_combined_first_130,R2_glms_combined_last_130,df_ready_dataset,df_outcome_covariates)

#Final multivariate regression analysis 

#Results table
results <- data.frame(Predictor=rownames(cols_predictors),numberofsig=NA,sig=NA)
results$numberofsig<-results_discovery$Discovery_significance+results_replication$Replication_significance
results$sig_in_both_datasets<-ifelse(results$numberofsig==2,1,0)
sigvariables<-results$Predictor[which(results$sig_in_both_datasets==1)]
#sigvariables will be used in volcano plotting. 

write.csv(results,file="results.csv",row.names = FALSE)

#Final multivariate analysis

outcome_and_covariate<-c("PE","sex_f31_0_0","age_at_recruitment_f21022_0_0")

final_multivariate_variables<-c(outcome_and_covariate,sigvariables)

df<-df[,grep(paste(final_multivariate_variables, collapse="|"),names(df))]

df<- df[complete.cases(df), ]
final_sample_size<-as.data.frame(nrow(df))
write.csv(final_sample_size,file="final_model_sample_size.csv")

final_multivariatemodel<-summary(glm(PE ~., data = df, family = binomial))
coeffinal_mul<-final_multivariatemodel$coefficients
final_multivariateresults <- data.frame(Predictor=rownames(coeffinal_mul),Final_Pval=NA,Final_Coefficient=NA,Final_SD=NA,Final_Oddsratios=NA)

#Evaluating the nagelkerke R2 of the final model. 
Final_R2<-nagelkerke(glm(PE ~., data = df, family = binomial))
R2_final=Final_R2$Pseudo.R.squared.for.model.vs.null[3,1]
Final_Nagelkerke_R2<-as.data.frame(R2_final)
write.csv(Final_Nagelkerke_R2,file="final_model_nagelkerke_pseudo_R2.csv",row.names = FALSE)

for(i in 1:nrow(final_multivariateresults)){
  final_multivariateresults[i,2] <- final_multivariatemodel$coefficients[i,4] ##line i,4.column
  final_multivariateresults[i,3] <- final_multivariatemodel$coefficients[i,1]
  final_multivariateresults[i,4] <- final_multivariatemodel$coefficients[i,2]
}
final_multivariateresults$Final_Oddsratios<- exp(final_multivariateresults$Final_Coefficient)
final_multivariateresults$Final_log10Pval<--log10(final_multivariateresults$Final_Pval)
final_multivariateresults$Final_significance<-ifelse(final_multivariateresults$Final_Pval<=(0.05),1,0)
final_multivariateresults$Final_ORlci <- exp(final_multivariateresults$Final_Coefficient + qnorm(c(0.025)) * final_multivariateresults$Final_SD)
final_multivariateresults$Final_ORhci <- exp(final_multivariateresults$Final_Coefficient + qnorm(c(0.975)) * final_multivariateresults$Final_SD)

write.csv(final_multivariateresults,file="final_multivariateresults.csv",row.names=FALSE)

remove(Final_R2,Final_Nagelkerke_R2,final_multivariate_variables,outcome_and_covariate,sigvariables,final_multivariatemodel,coeffinal_mul)

allresults <- data.frame(Predictor=rownames(cols_predictors),Discovery_Pval=NA,Discovery_OR=NA,Discovery_R2=NA,Replication_Pval=NA,Replication_OR=NA,Replication_R2=NA)
allresults$Discovery_Pval<-results_discovery$Discovery_Pval
allresults$Replication_Pval<-results_replication$Replication_Pval
allresults$Discovery_OR<-results_discovery$Discovery_Oddsratios
allresults$Replication_OR<-results_replication$Replication_Oddsratios
allresults$Discovery_R2<-results_discovery$Discovery_R2
allresults$Replication_R2<-results_replication$Replication_R2
allresults$sig_in_discovery_and_replication_datasets<-results$sig_in_both_datasets
allresults$sig_in_discovery_and_replication_datasets<-as.factor(allresults$sig_in_discovery_and_replication_datasets)
allresults$Discovery_significance<-ifelse(allresults$Discovery_Pval <=(0.05/nrow(cols_predictors)),1,0)
allresults$Replication_significance<-ifelse(allresults$Replication_Pval <=(0.05/nrow(cols_predictors)),1,0)
allresults$Combined_Pval<-results_combined$Combined_Pval
allresults$Combined_OR<-results_combined$Combined_Oddsratios
allresults$Combined_R2<-results_combined$Combined_R2

allresults_significant<-allresults[which(allresults$sig_in_discovery_and_replication_datasets=="1"),]

write.csv(allresults,file="allresults.csv",row.names = FALSE)
patterns=c("Yes","No")
for (pattern in patterns){allresults <-data.frame(lapply(allresults,function(x)str_replace_all(x,fixed(pattern),"")))}
allresults$Predictor<-lapply(allresults$Predictor,function(x)gsub("[_]"," ",x))

write_xlsx(allresults,"Supplementary_Material5_allresults.xlsx")

for (pattern in patterns){allresults_significant <-data.frame(lapply(allresults_significant,function(x)str_replace_all(x,fixed(pattern),"")))}
for (pattern in patterns){final_multivariateresults <-data.frame(lapply(final_multivariateresults,function(x)str_replace_all(x,fixed(pattern),"")))}
allresults_significant$Predictor<-lapply(allresults_significant$Predictor,function(x)gsub("[_]"," ",x))
final_multivariateresults$Predictor<-lapply(final_multivariateresults$Predictor,function(x)gsub("[_]"," ",x))

write_xlsx(allresults_significant,"Table2.xlsx")
write_xlsx(final_multivariateresults,"Table3.xlsx")

remove(allresults_significant,cols_predictors,df,final_multivariateresults,results,results_combined,results_discovery,results_replication,i,R2_final)

#Step 9: Plotting

#3 manhattan plots(Discovery, Replication and Combined(n=155,274)) and three volcano plots

#Upload the predictors' file, the information of groups(Biological samples, cognitive function ...) will be used in the Manhattan Plot

df_predictors <- data.table::fread("Predictor-table.csv",data.table = F,stringsAsFactors = F)

#Upload all results
allresults <- data.table::fread("allresults.csv",data.table = F,stringsAsFactors = F)

names_of_all_variables<-allresults$Predictor
#Turn significance into factor
allresults$sig_in_discovery_and_replication_datasets<-as.factor(allresults$sig_in_discovery_and_replication_datasets)

df_predictors$Predictor<-str_c(df_predictors$`Predictor Variable`,df_predictors$`f code`, sep="_f")
#Save the Supplementary table 4 
rownames(df_predictors)<-str_c(df_predictors$`Predictor Variable`,df_predictors$`f code`, sep="_f")
supptable4<- df_predictors[rownames(df_predictors)%in% names_of_all_variables,] 

write_xlsx(supptable4,"Supplementary_table_4_Predictors.xlsx")

df_predictors<-subset(df_predictors,select=c(Group,Predictor))

df<-merge(allresults,df_predictors,by="Predictor")
#Order by groups alphabetically
df <- df[with(df,order(Group)),]
#Variable position is assigned in ascending order, so that predictors in the same group are clustered together in the Manhattan Plot
#I make variable positions for the x-axis, and leave 3 position absent inbetween groups to highlight which variables belong to which group
df$Variable_position<-c(1:46,50:51,55:59,63:83,87:135,139:155,159,163,167:217,221:222,226:244,248:249,253:254,258:278,282,286:292)

final_multivariateresults <- data.table::fread("final_multivariateresults.csv",data.table = F,stringsAsFactors = F)
#Remove intercept 
final_multivariateresults <- final_multivariateresults[-which(final_multivariateresults$Predictor=="(Intercept)"),]
#Order variables in ascending order of their final_Pvalue
final_multivariateresults <- final_multivariateresults[with(final_multivariateresults,order(Final_Pval)),]
#Get the top three associated variants
top3variants_in_final_analysis <- final_multivariateresults$Predictor[1:3]
top1variant_in_final_analysis <- final_multivariateresults$Predictor[1]

patterns=c("Yes","No")
for (pattern in patterns){top3variants_in_final_analysis <-data.frame(lapply(top3variants_in_final_analysis,function(x)str_replace_all(x,fixed(pattern),"")))}
for (pattern in patterns){top1variant_in_final_analysis <-data.frame(lapply(top1variant_in_final_analysis,function(x)str_replace_all(x,fixed(pattern),"")))}

df_top<- df[which(df$Predictor %in% c(top3variants_in_final_analysis)),]
df_top1<- df[which(df$Predictor %in% c(top1variant_in_final_analysis)),]

df_top$Predictor<-lapply(df_top$Predictor,function(x)gsub("[_]"," ",x))
df_top1$Predictor<-lapply(df_top1$Predictor,function(x)gsub("[_]"," ",x))

names(df)<-lapply(names(df),function(x)gsub("[_]"," ",x))
names(allresults)<-lapply(names(allresults),function(x)gsub("[_]"," ",x))

#In discovery dataset

windowsFonts(Ariel=windowsFont("Ariel"))
#1
mp_discovery <- ggplot(df,aes(x=`Variable position`,y=-log10(`Discovery Pval`),color=Group)) + geom_point() + geom_hline(yintercept = -log10(.05/nrow(df)), col = c("black"), lwd = 0.5, linetype = c("dashed"))
mp_discovery<-mp_discovery + theme(axis.title.x = element_blank())
mp_discovery<-mp_discovery+ scale_x_continuous(breaks=NULL)
mp_discovery<-mp_discovery + theme(text=element_text(family="Ariel",size=12))

output_discovery = mp_discovery +
  geom_point(data = df_top, aes(x = Variable_position, y = -log10(Discovery_Pval)),  col = "gold", size = 5, shape = 18) +
  geom_text_repel(data = df_top, aes(label = Predictor, x= Variable_position, y =-log10(Discovery_Pval)),
                  nudge_y = 0.3, segment.size = 0.3, segment.colour = "grey", size = 3, col = "black")
output_discovery<-output_discovery + theme(text=element_text(family="Ariel",size=12))

jpeg(file="Manhattandiscovery.jpeg")
output_discovery
dev.off()

# for replication dataset
mp_replication <- ggplot(df,aes(x=`Variable position`,y=-log10(`Replication Pval`),color=Group)) + geom_point() + geom_hline(yintercept = -log10(.05/nrow(df)), col = c("black"), lwd = 0.5, linetype = c("dashed"))
mp_replication<-mp_replication + theme(axis.title.x = element_blank())
mp_replication<-mp_replication+ scale_x_continuous(breaks=NULL)
mp_replication<-mp_replication + theme(text=element_text(family="Ariel",size=12))

output_replication = mp_replication +
  geom_point(data = df_top, aes(x = Variable_position, y = -log10(Replication_Pval)),  col = "gold", size = 5, shape = 18) +
  geom_text_repel(data = df_top, aes(label = Predictor, x= Variable_position, y =-log10(Replication_Pval)),
                  nudge_y = 0.3, segment.size = 0.3, segment.colour = "grey", size = 3, col = "black")
output_replication<-output_replication + theme(text=element_text(family="Ariel",size=12))

jpeg(file="Manhattanreplication.jpeg")
output_replication
dev.off()

#Mp for Combined

mp_combined <- ggplot(df,aes(x=`Variable position`,y=-log10(`Combined Pval`),color=Group)) + geom_point() + geom_hline(yintercept = -log10(.05/nrow(df)), col = c("black"), lwd = 0.5, linetype = c("dashed"))
mp_combined<-mp_combined + theme(axis.title.x = element_blank())
mp_combined<-mp_combined+ scale_x_continuous(breaks=NULL)

output_Combined = mp_combined +
  geom_point(data = df_top, aes(x = Variable_position, y = -log10(Combined_Pval)),  col = "gold", size = 5, shape = 18) +
  geom_text_repel(data = df_top, aes(label = Predictor, x= Variable_position, y =-log10(Combined_Pval)),
                  nudge_y = 0.3, segment.size = 0.3, segment.colour = "grey", size = 3, col = "black")

jpeg(file="ManhattanCombined.jpeg")
output_Combined
dev.off()

#Volcano plot (discovery)

vp_discovery <- ggplot(allresults, aes(Discovery_OR, -log10(Discovery_Pval), 
                                       color=`sig in discovery and replication datasets`))
vp_discovery <- vp_discovery + geom_point(alpha=.5) + scale_x_log10() 
vp_discovery <- vp_discovery + scale_color_manual(breaks = c("1", "0"),
                                                  values=c("orange", "darkgrey"))

jpeg(file="VolcanoPlotDiscovery.jpeg")
vp_discovery
dev.off()

#Volcano plot (Replication)

vp_replication <- ggplot(allresults, aes(Replication_OR, -log10(Replication_Pval), 
                                         color=`sig in discovery and replication datasets`))
vp_replication <- vp_replication + geom_point(alpha=.5) + scale_x_log10() 
vp_replication <- vp_replication + scale_color_manual(breaks = c("1", "0"),
                                                      values=c("orange", "darkgrey"))

jpeg(file="VolcanoPlotReplication.jpeg")
vp_replication
dev.off()

#Volcano plot (combined)

vp_combined <- ggplot(allresults, aes(`Combined OR`, -log10(`Combined Pval`), 
                                      color=`sig in discovery and replication datasets`))
vp_combined <- vp_combined + geom_point(alpha=.5) + scale_x_log10() 
vp_combined <- vp_combined + scale_color_manual(breaks = c("1", "0"),
                                                values=c("orange", "darkgrey"))
vp_combined = vp_combined +
  geom_point(data = df_top1, aes(x = Combined_OR, y = -log10(Combined_Pval)),  col = "red", size = 7, shape = 18) +
  geom_text_repel(data = df_top1, aes(label = Predictor, x= Combined_OR, y =-log10(Combined_Pval)),
                  nudge_y = 0.3, segment.size = 0.3, segment.colour = "grey", size = 3, col = "black",force =5)

vp_combined<-vp_combined + theme(text=element_text(family="Ariel",size=12))

jpeg(file="VolcanoPlotCombined.jpeg")
vp_combined
dev.off()

remove(Ariel,df,vp_combined,vp_discovery,vp_replication,mp_replication,mp_combined,mp_discovery,output_discovery,output_Combined,output_replication,df_top,df_top1,top1variant_in_final_analysis,df_predictors,allresults,supptable4,top3variants_in_final_analysis,names_of_all_variables,find_rows_top,final_sample_size)