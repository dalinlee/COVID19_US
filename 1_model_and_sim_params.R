work_dir = getwd()
res_dir=paste(work_dir,"/results/",sep="")
dir.create(res_dir)

#### parameters across all models
inc_d=6      #### incubation period
lat_d = 2    #### latent period
ser_d=7.5    #### serial interaval

dis_d = 14
#### traffic parmaters for Wuhan area 
liw= 3546 # wuhan internationalinbound   
lwi=3633 # wuhan international outbound
lcw1 = 487310 #wuhan domestic inbound pre spring festiival
lcw2=810500 # wuhan domestic inbound post spring festival
lwc1=502013 # wuhan domestic outbound pre spring festival
lwc2=717226 #  wuhan domestic outbound post print festival 

N_us_out = 1.02E6
N_us_in = 1.05E6
#### date of flight from WUH to US
date_sfo=c(3,5,7,10,12,14,17,19,21,24,26,28,31,
           33,35,38,40,42,45,47,49,52,54)  +inc_d 
date_jfk = c(1,4,6,8,11,13,15,18,20,22,25,27,29,32,34,36,37,39,41,43,44,46,48,50,53)   + inc_d



flow=read.table("BaiduImmigration.csv",sep=",",header=T)
flow = flow[!is.na(flow$day),]
flow_sfo = flow[(flow$day + 6) %in% date_sfo,]
flow_jfk = flow[(flow$day + 6) %in% date_jfk,]
total_flow_sfo = sum(flow_sfo$index)
total_flow_jfk = sum(flow_jfk$index)

f_nostop_db1b=3338063
f_withstop_db1b=6508753
f_ratio=(f_nostop_db1b+f_withstop_db1b)/f_nostop_db1b   



