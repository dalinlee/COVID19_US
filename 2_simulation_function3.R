simulation_function <- function(sim, R0_red, R0, dis_d,
                                prop_z,ser_d,date_sim_end = "2020-04-01",
                                lat_d=2,inc_d=6,R0_red_108 = 0.5,
                                flow = flow,total_flow_sfo=total_flow_sfo,
                                total_flow_jfk = total_flow_jfk,f_ratio = f_ratio,
                                date_sfo=date_sfo,date_jfk=date_jfk,inf_count=T
                                ){
  # initializing intermediate parameters
  if(sim %% 10 ==0) print(sim)
  inits = para_init(ser_d=ser_d,lat_d = lat_d,
                    prop_z=prop_z,inc_d = inc_d,
                    R0=R0,date_sim_end = date_sim_end)
  res_wuhan = res_us=inits
  
  for (i in 1:inits$days_sim) {
    #print(i)
    #print(nrow(ind_aff))
    #### start simulating wuhan, till 01/23/2020
    if(i <= 54 + inc_d){ 
      res_wuhan = sim_Wuhan(i,inc_d,dis_d,err=inits$err,z_t = inits$z_t,
                           lcw1,lcw2=lcw2,lwc1,lwc2,
                           N=res_wuhan$N,S=res_wuhan$S,flow,total_flow_sfo,
                           total_flow_jfk, 
                           ind_aff=res_wuhan$ind_aff,n_inf_t=res_wuhan$n_inf_t,
                           us_imp=res_wuhan$us_imp,f_ratio,date_1st_in=res_wuhan$date_1st_in,
                           date_sfo,date_jfk)
      
     }
  ######## simulating the spreading in US.
    if(!is.null(res_wuhan$us_imp)){
       res_us = sim_us(i,err=inits$err,R0_red=R0_red,R0_red_108=R0_red_108,dis_d,
                      us_imp=res_wuhan$us_imp,ind_aff_us=res_us$ind_aff_us,
                      us_rec=res_us$us_rec,
                      N_us_out=N_us_out,N_us_in=N_us_in,
                      N_us=res_us$N_us,S_us=res_us$S_us,quara=res_us$quara,inf_count) 
     }
    
    
    
    
  }

  
  ######## output
  if(is.null(res_us$n_imp) | res_us$n_imp==0){
   n_imp=0
   date_1st_in=NA
   us_rec = rep(0,inits$days_sim - 60 +1)
   
  } else {
    n_imp = res_wuhan$n_imp_54
    date_1st_in = res_wuhan$date_1st_in
    us_rec = res_us$us_rec
    
  }
    
    
  us_seq_i= c(sim,R0,R0_red,prop_z,R0_red_108,
                n_imp,date_1st_in,
                us_rec)

  #us_seq = us_seq_i
  # print(us_seq)
  return(us_seq_i)
}




################### function for parameter initiation

para_init=function(ser_d,lat_d,prop_z,R0,inc_d,date_sim_end = "2020-04-01"){
  di=(ser_d - lat_d)*2   #### infections period
  dis_d = di + lat_d
  z=43*prop_z
  err = R0/di
  z_t= z/(31 + inc_d )
  
  # initialize population data
  N=1.9E7 # wuhan population
  S=N # wuhan sample
  N_us = 3.31E8 #us population
  S_us = N_us -1 # us sample
  
  
  # initialize number infected at begining of model
  n_inf_t = 1 # n infected at time t
  n_contag_t=1 # n contagious at time t 
  ind_aff <- data.frame("infected" = 1,"date_infected" = 0,"days_infected" = 1,"contagious" = 0)
  n_imp_44 =n_imp_46=n_imp_54=0
  n_us_98=n_us_113 = n_us_129 =NA
  us_imp={}
  ind_aff_us={}
  n_aff_us = 0
  n_inf_us_t = 0
  us_rec = {}
  days_sim = as.Date(date_sim_end) - as.Date("2019-11-24")
  date_1st_in = NA
  aff_wuhan_1st_in = NA
  quara={}
  n_inf_us_t=0
  res_us={}
  init_par=list(di,dis_d,z,err,z_t,N,S,N_us,S_us,
             n_inf_t,n_contag_t,ind_aff,n_imp_44,
             n_us_98,us_imp,ind_aff_us,n_aff_us,
             us_rec,days_sim,date_1st_in,aff_wuhan_1st_in,
             res_us,n_inf_us_t,quara)
  names(init_par) = c("di","dis_d","z","err","z_t","N","S","N_us","S_us",
                      "n_inf_t","n_contag_t","ind_aff","n_imp_44",
                      "n_us_98","us_imp","ind_aff_us","n_aff_us",
                      "us_rec","days_sim","date_1st_in","aff_wuhan_1st_in",
                      "res_us","n_inf_us_t","quara")
  return(init_par)
  
}



#########################################################
######### function for wuhan simulation
#########################################################
sim_Wuhan = function(i,inc_d,dis_d,err,z_t,
                     lcw1,lcw2,lwc1,lwc2,
                     N,S,flow ,total_flow_sfo,
                     total_flow_jfk , 
                     ind_aff,n_inf_t ,us_imp,
                     f_ratio,date_1st_in,date_sfo,date_jfk){
  
  
  
  ######## population immigration
  if(i<=47) {
    lcw = lcw1
    lwc = lwc1
  } else {
    lcw = lcw2
    lwc = lwc2
  }
  N_out=lwc+lwi
  N_in=liw + lcw
  
  
  prob_out = N_out/N
  n_aff = nrow(ind_aff)
  n_aff_out = rbinom(1,n_aff,prob_out)
  ids_aff =1:nrow(ind_aff)
  if (n_aff_out >0 & n_aff_out <= n_aff) {
    ids_aff_out = sample(ids_aff,n_aff_out)
    ind_aff = ind_aff[!ids_aff %in% ids_aff_out,]
    ids_aff=1:nrow(ind_aff)
  } 
  N=N - N_out + N_in
  S = S+ N_in - N_out * S/N - n_inf_t
  
  n_contag_t = sum(ind_aff$contag)
  
  
  ## infection
  #print(i)
  #print(paste("N of contagatious",n_contag))
  
  n_inf_exp= err *n_contag_t * S / N
  n_inf_t = rpois(1,n_inf_exp)
  n_inf_t = n_inf_t + ifelse(i<=31 + inc_d,rpois(1,z_t),0)
  #print(paste("N of infected at day i",n_inf_t))
  
  #n_inf_seq =c(n_inf_seq,n_inf_t)
  #id_inf_t = sample(1:N,n_inf_t)
  #n_inf_t=length(id_inf_t)
  
  if(n_inf_t >0) {
    ind_aff_t = cbind(rep(1,n_inf_t),rep(i,n_inf_t),rep(0,n_inf_t),rep(0,n_inf_t))
    colnames(ind_aff_t) = c("infected","date_infected","days_infected","contagious")
    ind_aff = rbind(ind_aff,ind_aff_t)  
  }
  #print(paste("total infected",nrow(ind_aff)))
  
  ##### export to US
  if(any(flow$day + inc_d == i)) flow_i = flow[flow$day == i -inc_d & !is.na(flow$day),"index"] else flow_i = 0
  if(i<=31 + inc_d) {
    n_fly_sfo = round(3239/13)
    n_fly_jfk = round(3209/13)
  } else {
    n_fly_sfo = round(4111 * 10/14* flow_i/total_flow_sfo )
    n_fly_jfk = round(3209 * 10/14* flow_i/total_flow_jfk )
    
  }
  
  #if(i<31 + 6) n_fly_sfo = round(3239/13) else n_fly_sfo = round(4111/14)
  #n_fly_jfk = round(3209/13)
  
  
  n_fly = ifelse(i %in% date_sfo, n_fly_sfo,ifelse(i %in% date_jfk,n_fly_jfk,0)) * f_ratio
  #n_fly = ifelse(i %in% date_sfo, n_fly_sfo,ifelse(i %in% date_jfk,n_fly_jfk,0)) 
  
  #print(n_fly)
  #ids_aff = 1:nrow(ind_aff)
  n_aff = nrow(ind_aff)
  
  
  if(n_fly>0){
    prob_fly = n_fly/N
    #print(i)
    #print(paste("totally",n_aff,"affected"))
    #print(paste("probably of fly",prob_fly))
    
    n_aff_fly = rbinom(1,n_aff,prob_fly)
    #print(paste("number of case fly",n_aff_fly))
    ids_aff = 1:nrow(ind_aff)
    if(n_aff_fly>0){
      ids_aff_fly= sample(ids_aff,n_aff_fly)
      #print(paste(length(ids_aff_fly),"cases move to us"))
      us_imp_i = ind_aff[ids_aff_fly,]
      if(is.null(us_imp)){
        date_1st_in = i
        aff_wuhan_1st_in = n_aff
      } 
      us_imp = rbind(us_imp,us_imp_i)
      us_imp=as.data.frame(us_imp)
      #ind_aff = ind_aff[ids_aff[!ids_aff %in% ids_us_aff],]
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

  if(i == 60 ) n_imp_54 = nrow(us_imp) else n_imp_54=NA
  Wuhan_res = list(N,S,us_imp=us_imp,ind_aff = ind_aff,n_inf_t = n_inf_t,n_imp_54,date_1st_in)
  names(Wuhan_res) = c("N","S","us_imp","ind_aff","n_inf_t","n_imp_54","date_1st_in")
  return(Wuhan_res)  
}





#########################################################
######### function for US simulation
#########################################################

sim_us = function(i,err,R0_red,R0_red_108,dis_d,
                  us_imp,ind_aff_us,us_rec=us_rec,
                  N_us_out,N_us_in,N_us,S_us,quara,inf_count=T){
  
  
  n_imp = nrow(us_imp)
  
  #### immigration in and out of US
  if(!is.null(ind_aff_us)){
    
    prob_out_us = N_us_out/N_us
    n_aff_us = ind_aff_us[,"N"]
    n_aff_us[n_aff_us<=0] = 0
    n_aff_us_out_exp = n_aff_us * prob_out_us
    
    n_aff_us_out = rpois(nrow(ind_aff_us),n_aff_us_out_exp)
    ind_aff_us[,"N"] = ind_aff_us[,"N"] - n_aff_us_out
    #print(paste(n_aff_us_out,"will be removed"))
    #id_aff_out = sample(1:nrow(ind_aff_us),n_aff_us_out)
    #ind_aff_us = ind_aff_us[-1 * id_aff_out,] 
  }
  if(is.null(ind_aff_us))n_inf_us_t = 0  else n_inf_us_t = sum(ind_aff_us[ind_aff_us[,"days_infected"]==1,"N"])
  N_us = N_us  - N_us_out + N_us_in
  S_us = S_us  - n_inf_us_t + N_us_in - N_us_out * S_us/N_us 
  
  
  
  ##### infection
  n_contag_imp = sum(us_imp[,"contagious"])
  if(!is.null(ind_aff_us)) n_contag_aff = sum(ind_aff_us[ind_aff_us[,"contagious"]==1,"N"],na.rm=T) else n_contag_aff  =0
  n_contag_us = n_contag_imp + n_contag_aff 
  if(n_contag_us<=0) n_contag_us = 0
  #print(paste("n contag in us",n_contag_us))
  #print(paste(n_contag_us,"are contagious"))
  n_inf_us_exp= err *n_contag_us * ifelse(i<=108,R0_red,R0_red_108*R0_red) * S_us/N_us 
  
  #print(paste(n_inf_us_exp,"are expected to be infected"))
  n_inf_us_t = rpois(1,n_inf_us_exp)
  #print(paste(n_inf_us_t,"are infected"))
  
  if(n_inf_us_t>0){
    #print(paste(nrow(ind_aff_us),"are previously infected in us"))
    
    #ind_aff_us_t = c(n_inf_us_t,1,i,0,0)
    #ind_aff_us_t = cbind(rep(1,n_inf_us_t),rep(i,n_inf_us_t),rep(0,n_inf_us_t),rep(0,n_inf_us_t))
    ind_aff_us_t = c(n_inf_us_t,1,i,0,0)
    names(ind_aff_us_t) = c("N","infected","date_infected","days_infected","contagious")
    ind_aff_us = rbind(ind_aff_us,ind_aff_us_t)    
    ind_aff_us = as.data.frame(ind_aff_us) 

    
  }
  
  
  
  ######## update the status in US cases  
  if(i > 54 + inc_d  ) {
    us_imp[,"days_infected"] = us_imp[,"days_infected"] + 1
    us_imp[us_imp[,"days_infected"]>=lat_d ,"contagious"] = 1
    us_imp[us_imp[,"days_infected"]>dis_d ,"contagious"] = 0
  }
  
  if(i >= 51 + inc_d){ us_imp[1,"contagious"] = 0 }
  if(i >= 55 + inc_d ) us_imp[1:min(nrow(us_imp),2),"contagious"] =0
  if(i >= 57 + inc_d ) us_imp[1:min(nrow(us_imp),5),"contagious"] =0
  if(i >= 62 + inc_d ) us_imp[1:min(nrow(us_imp),6),"contagious"] =0
  if(i >= 63 + inc_d ) us_imp[1:min(nrow(us_imp),7),"contagious"] =0
  if(i >= 64 + inc_d ) us_imp[1:min(nrow(us_imp),8),"contagious"] =0
  
  
  if(!is.null(ind_aff_us) ){
    ind_aff_us[,"days_infected"] = ind_aff_us[,"days_infected"] + 1
    ind_aff_us[ind_aff_us[,"days_infected"]>=lat_d ,"contagious"] = 1
    ind_aff_us[ind_aff_us[,"days_infected"]>dis_d ,"contagious"] = 0
    
    #if( i >=61 +inc_d ) ind_aff_us[1,"contagious"] =0
    if( i ==61 +inc_d ) {
      quara_t=ind_aff_us[1,]
      quara_t[1,c("N","contagious")] = c(1,0)
      quara = quara_t
      ind_aff_us[1,"N"] = ind_aff_us[1,"N"] -1
      #ind_aff_us = rbind(ind_aff_us,quara1)
    }
    
    if(nrow(ind_aff_us)==2 & i ==64 + inc_d ){
      quara_t=ind_aff_us[min(2,nrow(ind_aff_us)),]
      quara_t[1,c("N","contagious")] = c(1,0)
      quara = rbind(quara,quara_t)
      ind_aff_us[min(2,nrow(ind_aff_us)),"N"] = ind_aff_us[min(2,nrow(ind_aff_us)),"N"] -1
    } 
    if(i==89 + inc_d & !is.null(ind_aff_us)){
      quara_t=ind_aff_us[min(3,nrow(ind_aff_us)),]
      quara_t[1,c("N","contagious")] = c(2,0)
      quara = rbind(quara,quara_t)
      ind_aff_us[min(3,nrow(ind_aff_us)),"N"] = ind_aff_us[min(3,nrow(ind_aff_us)),"N"] -2            
      
    } 
    
    if(i==90 + inc_d & !is.null(ind_aff_us)) {
      
      quara_t=ind_aff_us[min(4,nrow(ind_aff_us)),]
      quara_t[1,c("N","contagious")] = c(3,0)
      quara = rbind(quara,quara_t)
      ind_aff_us[min(4,nrow(ind_aff_us)),"N"] = ind_aff_us[min(4,nrow(ind_aff_us)),"N"] -3           
    }
  }
  
  
  ########## count
  n_aff_us = ifelse(is.null(us_imp),0,sum(ind_aff_us[,"N"]))
  n_imp = ifelse(is.null(us_imp),0,nrow(us_imp))
  n_quara = ifelse(is.null(quara),0,sum(quara[,"N"]))
  
  n_us_t = n_aff_us+n_imp + n_quara    

  if (i >= 60 & !inf_count) us_rec = c(us_rec,n_us_t)
  if (i >= 60 & inf_count) us_rec = c(us_rec,n_contag_us)
  
  res_us = list(N_us,S_us,ind_aff_us,quara,n_imp,us_rec)
  names(res_us) = c("N_us","S_us","ind_aff_us","quara","n_imp","us_rec")
  return(res_us)
  
}
