#########  Estimate the scale of COVID-19 epidemic in US  #########
#########  As of 03/01/2020 using simulations ########## 
#########  Author: Dalin Li
#########  Date: 02/26/2020
#########  updated @ 03/16/2020

rm(list=ls(all=T))
work.dir=getwd()
setwd(work.dir)

source("1_model_and_sim_params.R")
source("2_simulation_function3.R")

library(parallel)


#### key parameters in simulation
n_sim = 1000 ## number of simulations
R0=2.3   #### R0, basic reproduction number
prop_z = 3 #### multiplier for z, the zoonotic force 
R0_red = 1 #### proportion of reduced R0 in U.S. for unidentified cases, before 03/11
R0_red_108 = 0.75 #### proportion of reduced R0 in U.S. for unidentified cases, after 03/11


us_seq <- mclapply(1:n_sim, simulation_function, R0_red, R0, dis_d,
                         prop_z,ser_d,date_sim_end = "2020-04-01",
                         lat_d=lat_d,inc_d=inc_d,R0_red_108 = R0_red_108,
                         flow = flow,total_flow_sfo=total_flow_sfo,
                         total_flow_jfk = total_flow_jfk,f_ratio = f_ratio,
                         date_sfo=date_sfo,date_jfk=date_jfk,inf_count=F,mc.cores = 44)
us_seq <- do.call(rbind, us_seq)
colnames(us_seq) = c("sim","R0","R0_red","prop_z","R0_red_108","n_imp_54","date_1st_in",
                           paste("n_contag_",60:129,sep=""))
      
us_seq =  as.data.frame(us_seq)
name_out = paste(res_dir,"us_seq_R0",R0,"_prop_z",prop_z,"_R0_red",
                       R0_red,"_R0_red_108",R0_red_108,"_03162020.txt",sep="")
write.table(us_seq,name_out,row.names=F,col.names=T,sep="\t")
 

