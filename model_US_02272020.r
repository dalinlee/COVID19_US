rm(list=ls(all=T))
par.dir<-"/mnt/share4/Dalin/projects/Corv/"
source("/mnt/share7/Guolab_share2/Dalin/projects/Simp_fun/simp.fun.r")
simp.create(par.dir,"")

R0 = 2
prop_z=1.5

inc_d=6
lat_d = 2
ser_d=7.5
di=(ser_d - lat_d)*2
dis_d = lat_d + di

liw= 3546
lwi=3633
lcw1 = 487310
lcw2=810500
lwc1=502013
lwc2=717226
date_sfo=c(3,5,7,10,12,14,17,19,21,24,26,28,31,
           33,35,38,40,42,45,47,49,52,54)  +6

date_jfk = c(1,4,6,8,11,13,15,18,20,22,25,27,29,32,34,36,39,41,43,46,48,50,53)   +6

n_sim = 1000
R0.c = c(2.1,2.2,2.3,2.4,2.5)
prop_z.c=c(1.5,2,2.5,3)
R0_red.c=c(1,0.75,0.5)

for(R0_red in R0_red.c){
  for (R0 in R0.c){
    for (prop_z in prop_z.c){
      print(R0)
      print(prop_z)
      print(R0_red)
      us_seq = {}  #### # of us cases imported in 01/23,total cases @ 01/23 and 03/01 (day 60,98 )
      for(sim in 1:n_sim){
        if(sim %% 10 ==0) print(sim)
        
        ########### initialize the intermediate parameters
        z=43*prop_z
        err = R0/di
        z_t= z/31
        n_inf=1
        n_contig =1
        n_diag=0
        N=1.9E7
        n_susp_t=N-1
        n_inf_t = 1
        n_contig_t=1
        n_diag_t=0
        
        n_inf_seq=1
        n_contig_seq = 1
        n_susp_seq = N-1
        lat_d=0
        T1=T2=T3=F
        ind_aff = matrix(c(1,0,1,1),nrow=1,ncol=4)
        ind_aff=as.data.frame(ind_aff)
        colnames(ind_aff) = c("infected","date_infected","days_infected","contagious")
        n_inf_seq=1
        us_imp={}
        ind_aff_us={}
        
        
        day_incub = 92    +6
        
        for (i in 1:day_incub){
          if(i<= 54    +6){
            id_aff =1:nrow(ind_aff)
            ######## population migration
            #N=nrow(ind_a)
            if(i<=41) {
              lcw = lcw1
              lwc = lwc1
            } else {
              lcw = lcw2
              lwc = lwc2
            }
            N_out=lwc+lwi
            N_in=liw + lcw
            
            ### population change
            id_out = sample(1:N,N_out,replace=F)
            n_out_case = sum(id_out %in% id_aff)
            if (n_out_case < nrow(ind_aff)& n_out_case>0) {
              id_aff_out = sample(id_aff,n_out_case)
              ind_aff = ind_aff[!id_aff %in% id_aff_out,]
              id_aff=1:nrow(ind_aff)
            }
            N=N - N_out + N_in
            
            ## infection
            #print(paste("day",i))
            #print(paste("contig individuals",n_contig_t))
            n_inf_exp= err *n_contig_t
            #print(paste("expected inf non-zoo",n_inf_exp))
            n_inf_t = rpois(1,n_inf_exp)
            #print(paste("observed inf non-zoo",n_inf_t))
            n_inf_t = n_inf_t + ifelse(i>6 & i<=31    +6,rpois(1,z_t),0)
            
            #print(paste("total inf,with zoo inf ",n_inf_t))
            n_inf_seq =c(n_inf_seq,n_inf_t)
            id_inf_t = sample(1:N,n_inf_t)
            #print(paste(sum(id_inf_t %in% id_aff),"previously infected"))
            id_inf_t = id_inf_t[!id_inf_t %in% id_aff]  
            n_inf_t=length(id_inf_t)
            #print(paste("final inf",n_inf_t))
            #print(paste("previous inf",nrow(ind_aff)))
            if(n_inf_t >0) {
              ind_aff_t = cbind(rep(1,n_inf_t),rep(i,n_inf_t),rep(0,n_inf_t),rep(0,n_inf_t))
              colnames(ind_aff_t) = c("infected","date_infected","days_infected","contagious")
              ind_aff = rbind(ind_aff,ind_aff_t)  
            }
           print(paste("total inf",nrow(ind_aff)))
            
            ##### export to US
            if(i<31 + 6) n_fly_sfo = round(3239/13) else n_fly_sfo = round(4111/14)
            n_fly_jfk = round(3209/13)
            n_fly = ifelse(i %in% date_sfo, n_fly_sfo,ifelse(i %in% date_jfk,n_fly_jfk,0))
            
            ids=1:N
            ids_aff = 1:nrow(ind_aff)
            if(n_fly>0){
              ids_us = sample(ids, n_fly)
              if(any(ids_us %in% ids_aff)){
                ids_us_aff= ids_us[ids_us %in% ids_aff]
                #print(paste(length(ids_us_aff),"cases move to us"))
                us_imp_i = ind_aff[ids_us_aff,]
                us_imp = rbind(us_imp,us_imp_i)
                us_imp=as.data.frame(us_imp)
                ind_aff = ind_aff[ids_aff[!ids_aff %in% ids_us_aff],]
              }
            }
            
            #### update
            if(nrow(ind_aff)>0){
              ind_aff[,"days_infected"] = ind_aff[,"days_infected"] + 1
              ind_aff[ind_aff[,"days_infected"]>=lat_d ,"contagious"] = 1
              ind_aff[ind_aff[,"days_infected"]>dis_d ,"contagious"] = 0
            }
            
            if(!is.null(us_imp) ){
              us_imp[,"days_infected"] = us_imp[,"days_infected"] + 1
              us_imp[us_imp[,"days_infected"]>=lat_d ,"contagious"] = 1
              us_imp[us_imp[,"days_infected"]>dis_d ,"contagious"] = 0
            }
            n_contig_t =sum(ind_aff[,"contagious"])
            
            
           # if(i == 52    +6){
           #    if(!is.null(us_imp)) T1 = T2 = T else T1=T2 = F
          #  }
            
           # if(i == 54    +6 & T1){
          #    if(nrow(us_imp)>5 ) T2 =T else T2= F
          #  }
            
          }   ###### finished the wuhan growth simulation
          
          
          if(!is.null(us_imp)>0 ){      
            #### calculating the totoal contagious individuals in US
            if(i >=54    +6 ) {
              us_imp[,"days_infected"] = us_imp[,"days_infected"] + 1
              us_imp[us_imp[,"days_infected"]>=lat_d ,"contagious"] = 1
              us_imp[us_imp[,"days_infected"]>dis_d ,"contagious"] = 0
              #us_imp[1:6,"contagious"] =0
            }
            
            if(i >= 51    +6){ us_imp[1,"contagious"] = 0 }
            if(i >=55    +6 ) us_imp[1:2,"contagious"] =0
            if(i >=57    +6 ) us_imp[1:min(nrow(us_imp),5),"contagious"] =0
            if(i >=62    +6 ) us_imp[1:min(nrow(us_imp),6),"contagious"] =0
            if(i >=63    +6 ) us_imp[1:min(nrow(us_imp),7),"contagious"] =0
            if(i >=64    +6 ) us_imp[1:min(nrow(us_imp),8),"contagious"] =0
            n_contig_imp = sum(us_imp[,"contagious"])
            
            if(!is.null(ind_aff_us) ) n_contig_dom = sum(ind_aff_us[,"contagious"]) else n_contig_dom = 0
              
            n_contig_us = n_contig_imp + n_contig_dom     
        
            #### infection
            n_inf_us_exp= err *n_contig_us * R0_red
            n_inf_us_t = rpois(1,n_inf_us_exp)
            #print(i)
            #print(paste("newly infected us cases",n_inf_us_t))
            if(n_inf_us_t>0){
              ind_aff_us_t = cbind(rep(1,n_inf_us_t),rep(i,n_inf_us_t),rep(0,n_inf_us_t),rep(0,n_inf_us_t))
              colnames(ind_aff_us_t) = c("infected","date_infected","days_infected","contagious")
              ind_aff_us = rbind(ind_aff_us,ind_aff_us_t)    
              ind_aff_us = as.data.frame(ind_aff_us)     
            }
            
            if(!is.null(ind_aff_us) ){
              ind_aff_us[,"days_infected"] = ind_aff_us[,"days_infected"] + 1
              ind_aff_us[ind_aff_us[,"days_infected"]>=lat_d ,"contagious"] = 1
              ind_aff_us[ind_aff_us[,"days_infected"]>dis_d ,"contagious"] = 0
              if( i >=61    +6 ) ind_aff_us[1,"contagious"] =0
              if(nrow(ind_aff_us)>2 & i >=64    +6 ) ind_aff_us[2,"contagious"] =0 
  
            }
            
            if(i == 54    +6) {
              n_imp = nrow(us_imp)
              n_us_60 = ifelse(!is.null(ind_aff_us),nrow(ind_aff_us) + n_imp,n_imp)
            }
            
            if( i == 92    +6 ){
              n_us_98 = ifelse(!is.null(ind_aff_us),nrow(ind_aff_us) + n_imp,n_imp)
            }
            
          }
        }
        if(is.null(us_imp)) n_imp = 0 else n_imp = nrow(us_imp)
        if(T1 & T2 ) us_seq_i= c(sim,R0,R0_red,prop_z,n_imp,n_us_60,n_us_98) else us_seq_i=c(sim,R0,R0_red,prop_z,n_imp,NA,NA)
        us_seq = rbind(us_seq,us_seq_i)
      }
      colnames(us_seq) = c("sim","R0","R0_red","prop_z","n_imp","n_us_60","n_us_98")
      simp.write(us_seq,paste(res.dir,"us_seq_R0",R0,"_prop_z",prop_z,"_R0_red",R0_red,".txt",sep=""))
    }
  }
} 
      
        






