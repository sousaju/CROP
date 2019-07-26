###################################
# Installing and loading packages #
###################################
pckgs=c("gplots","tidyverse","igraph","ggraph","robCompositions","compositions","gridExtra")
for (i in 1:length(pckgs))
{
  if (system.file(package = pckgs[i])=="") install.packages(pckgs[i]) 
}

if (system.file(package = "corrr")=="") {
  install.packages("devtools")  
  devtools::install_github("drsimonj/corrr")}

require(ggraph)
require(gplots)
require(tidyverse)
require(corrr)
require(igraph)
require(robCompositions)
require(compositions)
require(gridExtra)

#########################################
# Function for multiplicities detection #
#########################################
mltpl_detect=function (data,ccth,rtw,rt){
  set.seed(1)
  cor_tbl <- data %>% correlate() %>% stretch() %>% 
      filter(r >= ccth)
  
  rtx=numeric()
  rty=numeric()
  k=1
  row_ind=numeric()
  for (i in 1:length(cor_tbl$x))
  {
    rtx[i]=rt[1,which(names(rt) == cor_tbl$x[i])]
    rty[i]=rt[1,which(names(rt) == cor_tbl$y[i])]
    if (abs(rtx[i]-rty[i])<= rtw) {row_ind[k]=i; k=k+1}                                          
  }
  
  cor_tbl <- cor_tbl[row_ind,]
  
  return(cor_tbl)
}


####################################
# Function for correlation network #
####################################
mltpl_depict=function (det,ccth,rtw,stretched=NULL){  

  graph_corr <- det %>%  
    graph_from_data_frame(directed = FALSE)   #det is cor_tbl from previous function
if (class(stretched)!="NULL") { 
  clrs_str=ifelse(stretched=="no",rgb(255,255,255,maxColorValue=255),rgb(255,223,0,maxColorValue=255))
  cor_mapa = ggraph(graph_corr) +
    geom_edge_link(aes(edge_alpha = r, edge_width = r, color = r)) +
    scale_edge_colour_gradientn(limits = c(ccth, 1), colors = c(rgb(252,146,114,maxColorValue=255), "firebrick2")) +
    geom_node_point(aes(fill = clrs_str), size = 4,shape=21) + 
    scale_fill_identity("stretched \ncluster", labels = c("yes","no"), breaks = c(rgb(255,223,0,maxColorValue=255),rgb(255,255,255,maxColorValue=255)), guide = "legend") +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_graph() +
    labs(title = paste("CROP - Correlation-based Removal of Multiplicities
Correlation network with correlation coefficient above ",ccth," and retention time window +- ",rtw, sep=""))
} else {
  cor_mapa = ggraph(graph_corr) +
    geom_edge_link(aes(edge_alpha = r, edge_width = r, color = r)) +
    scale_edge_colour_gradientn(limits = c(ccth, 1), colors = c(rgb(252,146,114,maxColorValue=255), "firebrick2")) +
    geom_node_point(fill = "white", size = 4,shape=21) + 
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_graph() +
    labs(title = paste("CROP - Correlation-based Removal of Multiplicities
Correlation network with correlation coefficient above ",ccth," and retention time window +- ",rtw, sep=""))
}
  return(cor_mapa)
}

#######################################
# Function for multiplicities removal #
#######################################
mltpl_remove=function (data,det,name="",ccth,rtw,mcs,rt,maxrtw=NULL,funit,rtunit){
  # Selection of clusters
  cor_tbl2 <- det
  k=1
  max_clust_size = mcs
  zp <- tibble(c(rep("",max_clust_size)))
  nms_tot=character()
  index=1
  
  while(length(which(!(cor_tbl2$x %in% nms_tot))) > 0)
  {
    q = which(!(cor_tbl2$x %in% nms_tot))[1]
    mltpl <- filter(cor_tbl2, (x == cor_tbl2$x[q]) | (y == cor_tbl2$x[q]) | (x == cor_tbl2$y[q]) | (y == cor_tbl2$y[q]))
    nms =character()
    nms[1]=cor_tbl2$x[q]
    nms[2]=cor_tbl2$y[q]
    k=2
    while ((length(which(!(mltpl$x %in% nms))) > 0) | (length(which(!(mltpl$y %in% nms))) > 0))
    {
      for (i in which(mltpl$x %in% nms))
      {
        if(!(mltpl$y[i] %in% nms)) {
          mltpl <- rbind(mltpl,filter(cor_tbl2, ((x == mltpl$y[i]) & !(y %in% nms)) | ((y == mltpl$y[i]) & !(x %in% nms))))
          k=k+1
          nms[k]=mltpl$y[i]}
      }
      for (i in which(mltpl$y %in% nms))
      {
        if(!(mltpl$x[i] %in% nms)) {
          mltpl <- rbind(mltpl,filter(cor_tbl2, ((x == mltpl$x[i]) & !(y %in% nms)) | ((y == mltpl$x[i]) & !(x %in% nms))))
          k=k+1
          nms[k]=mltpl$y[i]}
      }
    }
    nms_tot = c(nms_tot,unique(nms))
    if (max_clust_size < length(unique(nms))) stop ("creating too big cluster(s)! \n      consider changing your parameteres ccth and/or rtw \n     not recommended but possible to override by setting mcs >= 200 \n  are you setting rtw in correct units ['min']?")
    cluster<-c(unique(nms),rep("",max_clust_size-length(unique(nms))))
    zp <- mutate(zp,cluster)
    zp <- zp %>% rename(!!paste('cluster',index,sep="") := cluster)
    index=index+1
    rm(mltpl)
  }
  
  zp <- zp[,-1]
  rm(nms)
  rm(cluster)
  rm(index)
  rm(k)
  rm(q)
  
  # Detection of stretched clusters
  if (class(maxrtw) != "NULL") {
    had_pos=numeric()
    had_con=numeric()
    for (i in 1:dim(zp)[2])
    {
      j=1
      while ((j <= dim(zp)[1]) & (as.character(zp[j,i])!=""))
      {
        had_con[j]=which(colnames(data)==as.character(zp[j,i])) 
        j=j+1
      }
      minrtp=had_con[which(rt[had_con]==min(rt[had_con]))]
      maxrtp=had_con[which(rt[had_con]==max(rt[had_con]))]
      if(rt[maxrtp[1]]-rt[minrtp[1]]>maxrtw) had_pos=c(had_pos,had_con)
      had_con=numeric()
    }
    
    mapa_pos=sort(which(colnames(data) %in% unique(c(det$x,det$y))))
    str_groups=rep("no",length(mapa_pos))
      
    for (k in 1:length(had_pos))
      str_groups[which(mapa_pos == had_pos[k])]="yes"
  }
  
  # Choice of representantion of clusters 
  i=1
  j=1
  k=1
  
  row_ind=colnames(data)
  data_mod=rbind(row_ind,rt,data)
  clustn_short=numeric()
  clustn_long=numeric()
  for (i in 1:dim(zp)[2])
  {
    mean = numeric()
    clmn_del = numeric()
    for (j in 1:(dim(unique(zp[,i]))[1]-1))
    {
      mean[j]=mean(as.numeric(data_mod[-c(1,2),which(colnames(data_mod) == as.character(zp[j,i]))]))
      names(mean)[j]=colnames(data_mod)[which(colnames(data_mod) == as.character(zp[j,i]))]
      clmn_del[j] = which(colnames(data_mod) == as.character(zp[j,i]))
    }

    colnames(data_mod)[clmn_del]
    clmn_keep = clmn_del[which.max(mean)]
    clmn_del = clmn_del[-which.max(mean)]
    nms_mltpl = colnames(data_mod)[clmn_keep]
    for (k in 1:length(clmn_del))
      nms_mltpl = paste(nms_mltpl,colnames(data_mod)[clmn_del[k]], sep=" + ")
    colnames(data_mod)[clmn_keep] = paste(colnames(data_mod)[clmn_keep],"*",sep="")
    data_mod[1,clmn_keep]=nms_mltpl
    clustn_short[i]=colnames(data_mod)[clmn_keep]
    clustn_long[i]=data_mod[1,clmn_keep]
    data_mod = data_mod[,-c(clmn_del)]
  }
  if (rtunit=="min") rownames(data_mod)=c("All Cluster Members","RT[min]",rownames(data)) else rownames(data_mod)=c("All Cluster Members","RT[s]",rownames(data))
  
  data_mod2=data_mod
  data_mod2[,ncol(data_mod)+1]=" "
  colnames(data_mod2)=c("Cluster Representation*",colnames(data_mod))
  write.table(data_mod2,paste(name,"_CROPped_ccth",ccth,"_rtw+-",rtw,"_data_without_multiplicities.csv",sep=""),sep=";")  
  
  otp_table=t(data_mod[-c(1,2),])
  otp_table=cbind(rownames(otp_table),otp_table)
  colnames(otp_table)[1]=funit
  rownames(otp_table)=NULL
  write.table(otp_table,paste(name,"_CROPped_ccth",ccth,"_rtw+-",rtw,"_output_table.csv",sep=""),sep=";",row.names=FALSE) 
  
  ls_clust=as.matrix(cbind(clustn_short,clustn_long,rep(" ",length(clustn_short))))
  rownames(ls_clust)=ls_clust[,1]
  ls_clust=ls_clust[,-1]
  colnames(ls_clust)=c("Cluster Representation","All Cluster Members")
  write.table(ls_clust,paste(name,"_CROPped_ccth",ccth,"_rtw+-",rtw,"_list_of_clusters.csv",sep=""),sep=";")
  
  if (class(maxrtw) != "NULL") return(str_groups) #vector with stretched clusters coding
}

#########################################
# CROP - Function for the whole process #
#########################################
CROP=function(data, name="", ccth = 0.75, rtw = 0.02, mcs=100, maxrtw = NULL, rtunit="min", funit="MW") {
  dataload = read.csv2(data,header=TRUE,dec=".")
  
  if (rtunit!="min" & rtunit!="s") stop("incorrect setting of rtunit   ['min' or 's']")
  
  if (dim(dataload)[2]==1) dataload = read.csv2(data,header=TRUE,dec=".",sep=",")
    
  if (dim(dataload)[2]==1) {
    stop ("only one column detected in the data \n   1: are you using a correct delimiter?   [';' or ','] \n   2: do you have first two columns with feature names \n      and retention time in the data?")
  } else if(dim(dataload)[2]<4) {
    stop("cannot compute correlations on 1 observation!")
  } else if(dim(dataload)[2]<7) {
    warning(paste("computing correlations only on ",dim(dataload)[2]-2," observations \n     consider using bigger sample size",sep=""))
  }
  
  if (dim(dataload)[1]==1) {
    stop("cannot compute correlations on 1 variable!")
  } else if ((dim(dataload)[1]<101) & (dim(dataload)[1]<dim(dataload)[2])) {
    warning(paste("computing correlations only on ",dim(dataload)[1]," variables \n     consider using more features \n  are you sure you have features in rows and samples in columns?",sep=""))
  } else if (dim(dataload)[1]<101) {
    warning(paste("computing correlations only on ",dim(dataload)[1]," variables \n     consider using more features",sep=""))
  } else if (dim(dataload)[1]<dim(dataload)[2]) {
    warning("are you sure you have features in rows and samples in columns?")
  }
  
  if (class(maxrtw) == "NULL") warning("you may have stretched clusters phenomenon in your data \n     which could not be highlighted in the correlation network \n  consider setting maxrtw for this purpose")
  
  rt_vec=dataload[,2]    
  rownames(dataload)=dataload[,1]
  dataload=dataload[,-c(1,2)] 
  dataload=t(dataload)
  rownames(dataload)=substring(rownames(dataload), 2)
  names(rt_vec)=colnames(dataload)
  rt_vec=as.data.frame(t(rt_vec))
  dataload=as.data.frame(dataload)
  
  if (min(dataload)<0) warning("are you using non normalized data? \n     negative values detected")
  
  if ((rtw > 0.1) & (rtunit=="min")) {
    warning("are you sure rtw is set in correct units? \n     current setting: rtunit = 'min' \n  possibly too large retention time window")
  } else if ((rtw < 0.6)& (rtunit=="s")) {
    warning("are you sure rtw is set in correct units? \n     current setting: rtunit = 's' \n  possibly too small retention time window")
  } else if ((max(rt_vec)>35) & (rtunit=="min")) {
    warning("are you sure retention time units are set correctly? \n     current setting: rtunit = min")
  } else if ((max(rt_vec)<180) & (rtunit=="s")) {
    warning("are you sure retention time units are set correctly? \n     current setting: rtunit = s")
  }
  
  # detection, depicting and removal of multiplicities
  data_pro_mapu=mltpl_detect(data=dataload, ccth=ccth, rtw=rtw, rt=rt_vec)
  
  if (class(maxrtw)!="NULL") {
    sc_identif=mltpl_remove(data=dataload, det=data_pro_mapu, name=name, ccth=ccth, rtw=rtw, mcs=mcs, rt=rt_vec, maxrtw=maxrtw,rtunit=rtunit,funit=funit)
    cairo_pdf(paste(name=name,"_CROPped_ccth",ccth,"_rtw+-",rtw,"_correlation_network.pdf",sep=""),w=25,h=25)
    print(mltpl_depict(det=data_pro_mapu, ccth=ccth, rtw=rtw, stretched=sc_identif))
    dev.off()
  } else {
    mltpl_remove(data=dataload, det=data_pro_mapu, name=name, ccth=ccth, rtw=rtw, mcs=mcs, rt=rt_vec, rtunit=rtunit,funit=funit)
    cairo_pdf(paste(name=name,"_CROPped_ccth",ccth,"_rtw+-",rtw,"_correlation_network.pdf",sep=""),w=25,h=25)
    print(mltpl_depict(det=data_pro_mapu, ccth=ccth, rtw=rtw))
    dev.off()
  }
}

###########################################################################################################################################
###########################################################################################################################################
##################
# CROP IN ACTION #
##################
CROP(data="example_data.csv", name="project_name", ccth=0.75, rtw=0.02, maxrtw=0.04)

###########################################################################################################################################
###########################################################################################################################################
#############################
# DESCRIPTION OF PARAMETERS #
#############################

###
# ccth ... threshold for correlation coefficient values 
#          default: ccth = 0.75
###
# funit ... units in which your features are measured ("MW", "m/z" or anything else)
#           only for column headers in the outputs
#           default: funit = "MW"
###
# maxrtw... maximal allowed RT window to color stretched clusters phenomenon
#                    recommended to set as max 2*rtw
#                    default: maxrtw = NULL
###
# mcs ... maximal allowed cluster size
#         should not be set high; if default is not enough, rather consider setting ccth bigger and/or rtw smaller than changing mcs
#         default: mcs = 100
###
# data ... filepath to an input csv table of your MS data with:
#          names of features in the first column
#          retention time of features in minutes in second column
#          samples in the rest of columns (including QCs)
#          column header as names of samples
###
# name ... a note which will be in names of all files
###
# rtunit ... units of retention time you ar using in your data ("min" or "s")
#            default: rtunit = "min"
###
# rtw ... retention time window where +-rtw will be considered
#         default: rtw = 0.02



