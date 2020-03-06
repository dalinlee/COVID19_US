#########  Estimate the scale of COVID-19 epidemic in US  #########
#########  As of 03/01/2020 using simulations ########## 
#########  Author: Dalin Li
#########  Date: 02/26/2020
#########  Clean up @ 03/05/2020


rm(list=ls(all=T))

res_dir = getwd()

#### parameters across all models
inc_d=6      #### incubation period
lat_d = 2    #### latent period
ser_d=7.5    #### serial interaval

di=(ser_d - lat_d)*2   #### infections period
dis_d = lat_d + di     #### disease duration(till non-infectious)


#### traffic parmaters for Wuhan area 
liw= 3546    
lwi=3633
lcw1 = 487310
lcw2=810500
lwc1=502013
lwc2=717226


#### date of flight from WUH to US
date_sfo=c(3,5,7,10,12,14,17,19,21,24,26,28,31,
         33,35,38,40,42,45,47,49,52,54)  +inc_d 
date_jfk = c(1,4,6,8,11,13,15,18,20,22,25,27,29,32,34,36,39,41,43,46,48,50,53)   + inc_d




#### key parameters in simulation
n_sim = 1000
R0_c = c(2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7)    #### R0, basic reproduction number
prop_z_c=c(1.5,2,2.5)                #### multiplier for z, the zoonotic force from the Huanai Seafood Market
R0_red_c=c(1,0.75)                 #### reduced R0 in U.S. for unidentified cases


### Now the simuation starts
for(R0_red in R0_red_c){
  for (R0 in R0_c){
    for (prop_z in prop_z_c){
      print(R0)
      print(prop_z)
      print(R0_red)
      us_seq = {}
      
      for(sim in 1:n_sim){

        if(sim %% 10 ==0) print(sim)
        z=43*prop_z
        err = R0/di
        z_t= z/31
        n_inf=1
        n_contig =1
        n_diag=0
  
        n_inf_t = 1
        n_contig_t=1
        
        N=19000000
        N_us=331000000
        S=N-1
  
        n_inf_seq=1
        n_contig_seq = 1
        n_susp_seq = N-1
        
        n_us_60=NA
        n_us_98=NA
       
        ind_aff = matrix(c(1,0,1,1),nrow=1,ncol=4)
        ind_aff=as.data.frame(ind_aff)
        colnames(ind_aff) = c("infected","date_infected","days_infected","contagious")


        n_inf_us=0
        us_imp={}
        ind_aff_us={}
        day_sim = 92 + inc_d
        
        

        for (i in 1:day_sim){
          #print(i)
          ##### simulate the spread in Wuhan area based on the SEIR model
          if(i<= 54 + inc_d){    ##### limit it to 01/23/2020 
          id_aff =1:nrow(ind_aff)
          if(i<=41 ) {
           lcw = lcw1
           lwc = lwc1
          } else {
           lcw = lcw2  
           lwc = lwc2
          } 
  
          N_out=lwc+lwi
          N_in=liw + lcw
          
          ### population change
          #ids=1:nrow(ind_a)
          S = S- n_inf_t + N_in - S/N*N_out
          N=N - N_out + N_in
          
          #### update the status of infected/contagious in Wuhan area
          if(nrow(ind_aff)>0){
            ind_aff[,"days_infected"] = ind_aff[,"days_infected"] + 1
            ind_aff[ind_aff[,"days_infected"]>=lat_d ,"contagious"] = 1
            ind_aff[ind_aff[,"days_infected"]>dis_d ,"contagious"] = 0
            if(i<6+lat_d)ind_aff[1 ,"contagious"] = 1
          }
          
          n_contig =sum(ind_aff[,"contagious"])
          
          n_inf_exp= err *n_contig * S/N
          n_inf_t = rpois(1,n_inf_exp)
  
          n_z = rpois(1,z_t)
          n_inf_t = n_inf_t + ifelse(i<=31  + inc_d,n_z,0)
          n_inf_seq =c(n_inf_seq,n_inf_t)
          id_inf_t = sample(1:N,n_inf_t)
  
          id_inf_t = id_inf_t[!id_inf_t %in% id_aff]
          n_inf_t=length(id_inf_t)
  
          if(n_inf_t >0) {
            ind_aff_t = cbind(rep(1,n_inf_t),rep(i,n_inf_t),rep(0,n_inf_t),rep(0,n_inf_t))
            colnames(ind_aff_t) = c("infected","date_infected","days_infected","contagious")
            ind_aff = rbind(ind_aff,ind_aff_t)
          }

          ##### export to US
          if(i<31 + inc_d) n_fly_sfo = round(3239/13) else n_fly_sfo = round(4111/14)
          n_fly_jfk = round(3209/13)
          n_fly = ifelse(i %in% date_sfo, n_fly_sfo,ifelse(i %in% date_jfk,n_fly_jfk,0))

  
          ids=1:N
          ids_aff = 1:nrow(ind_aff)
          if(n_fly>0){
            ids_us = sample(ids, n_fly)
   
            if(any(ids_us %in% ids_aff)){
              ids_us_aff= ids_us[ids_us %in% ids_aff]
              us_imp_i = ind_aff[ids_us_aff,]
              us_imp = rbind(us_imp,us_imp_i)
              us_imp=as.data.frame(us_imp)
    
              ind_aff = ind_aff[ids_aff[!ids_aff %in% ids_us_aff],]
    
              } 
            }

          
  
          }  ######## done with the simulation of disease spread in Wuhan 
   

          #### simulating disease spread in U.S.
          #### update the status of imported cases in U.S.    
          
          if(!is.null(us_imp) ){
            #print(i)
            #print(us_imp)
            us_imp[,"days_infected"] = us_imp[,"days_infected"] + 1
            us_imp[us_imp[,"days_infected"]>=lat_d ,"contagious"] = 1
            us_imp[us_imp[,"days_infected"]>dis_d ,"contagious"] = 0
            
            n_imp = nrow(us_imp)
            #### update the state of affeced US sample
            if(!is.null(ind_aff_us) ){
              ind_aff_us[,"days_infected"] = ind_aff_us[,"days_infected"] + 1
              ind_aff_us[ind_aff_us[,"days_infected"]>=lat_d ,"contagious"] = 1
              ind_aff_us[ind_aff_us[,"days_infected"]>dis_d ,"contagious"] = 0
              
              if( i >=61    + inc_d ) ind_aff_us[1,"contagious"] =0  
              if(nrow(ind_aff_us)>2 & i >=64    + inc_d ) ind_aff_us[2,"contagious"] =0   
            }  
            
            #### include individual events
            
            if(i >= 52    + inc_d){ us_imp[1,"contagious"] = 0 }  
            if(i >=55    + inc_d ) us_imp[1:min(nrow(us_imp),2),"contagious"] =0  
        
            if(i >=57    + inc_d ) us_imp[1:min(nrow(us_imp),5),"contagious"] =0  
            if(i >=62    + inc_d ) us_imp[1:min(nrow(us_imp),6),"contagious"] =0  
            if(i >=63    + inc_d ) us_imp[1:min(nrow(us_imp),7),"contagious"] =0  
            if(i >=64    + inc_d ) us_imp[1:min(nrow(us_imp),8),"contagious"] =0  
            
            if(i>=61 + inc_d & !is.null(ind_aff_us)) ind_aff_us[1,"contagious"] = 0
            if(i>=64 + inc_d & !is.null(ind_aff_us)) ind_aff_us[1:min(nrow(ind_aff_us),2),"contagious"] = 0
            if(i>=89 + inc_d & !is.null(ind_aff_us)) ind_aff_us[1:min(nrow(ind_aff_us),4),"contagious"] = 0
            if(i>=90 + inc_d & !is.null(ind_aff_us)) ind_aff_us[1:min(nrow(ind_aff_us),7),"contagious"] = 0
            
            n_contig_imp = sum(us_imp[,"contagious"])
            if(!is.null(ind_aff_us) ) n_contig_dom = sum(ind_aff_us[,"contagious"]) else n_contig_dom = 0
            n_contig_us = n_contig_imp + n_contig_dom
    
            ### susceptible  
            if(!is.null(ind_aff_us))S_us = N_us - n_inf_us else S_us = N_us
            
            #### infection
            n_inf_us_exp= err *n_contig_us * R0_red * S_us/N_us
            n_inf_us_t = rpois(1,n_inf_us_exp)
           if(n_inf_us_t>0){
              ind_aff_us_t = cbind(rep(1,n_inf_us_t),rep(i,n_inf_us_t),rep(0,n_inf_us_t),rep(0,n_inf_us_t))
              colnames(ind_aff_us_t) = c("infected","date_infected","days_infected","contagious")
              ind_aff_us = rbind(ind_aff_us,ind_aff_us_t)
              ind_aff_us = as.data.frame(ind_aff_us)
              n_inf_us =  nrow(ind_aff_us)
            } 
    
             ####### count  
             if(i == 54    + inc_d) {
               n_imp = nrow(us_imp)
               n_us_60 = ifelse(!is.null(ind_aff_us), n_inf_us + n_imp,n_imp)
             } 
             if( i == 92    + inc_d ){
               n_us_98 = ifelse(!is.null(ind_aff_us),n_inf_us + n_imp,n_imp)
            }

          } else n_imp = 0
          
        } 
        if(!is.null(ind_aff_us) ) us_seq_i= c(sim,R0,R0_red,prop_z,n_imp,n_us_60,n_us_98) else us_seq_i=c(sim,R0,R0_red,prop_z,n_imp,NA,NA)
        #us_seq_i
        us_seq = rbind(us_seq,us_seq_i) 
      }
      
      
    colnames(us_seq) = c("sim","R0","R0_red","prop_z","n_imp","n_us_60","n_us_98")
    write.table(us_seq,paste(res_dir,"us_seq_R0",R0,"_prop_z",prop_z,"_R0_red",R0_red,".txt",sep=""),row.names=F,col.names=T,sep="\t")

    }
  }
}







