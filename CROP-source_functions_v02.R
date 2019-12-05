###################################
# Installing and loading packages #
###################################
pckgs=c("gplots","tidyverse","igraph","ggraph","robCompositions","compositions","gridExtra","stringi")
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
require(stringi)

#########################################
# Function for multiplicities detection #
#########################################
mltpl_detect=function (data,ccth,rtw,rt){
  set.seed(1)
  cor_tbl <- suppressMessages(data %>% correlate() %>% stretch() %>% 
      filter(r >= ccth))
  
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
mltpl_depict=function (det,ccth,rtw,stretched=NULL,rtunit,mztab){  
  
  if (mztab==TRUE){
    for (i in 1:nrow(det)) {
      det$x[i]=strsplit(det$x,"|",fixed=TRUE)[[i]][2]
      det$y[i]=strsplit(det$y,"|",fixed=TRUE)[[i]][2]
    }}
  
  graph_corr <- det %>%  
    graph_from_data_frame(directed = FALSE)   #det is cor_tbl from previous function
if (class(stretched)!="NULL") { 
  clrs_str=ifelse(stretched=="no",rgb(255,255,255,maxColorValue=255),rgb(255,223,0,maxColorValue=255))
  cor_mapa = ggraph(graph_corr, layout='nicely') +
    geom_edge_link(aes(edge_alpha = r, edge_width = r, color = r)) +
    scale_edge_colour_gradientn(limits = c(ccth, 1), colors = c(rgb(252,146,114,maxColorValue=255), "firebrick2")) +
    geom_node_point(aes(fill = clrs_str), size = 4,shape=21) + 
    scale_fill_identity("stretched \ncluster", labels = c("yes","no"), breaks = c(rgb(255,223,0,maxColorValue=255),rgb(255,255,255,maxColorValue=255)), guide = "legend") +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_graph() +
    labs(title = paste("CROP - Correlation-based Removal of Multiplicities
Correlation network with correlation coefficient above ",ccth," and retention time window +- ",rtw, rtunit, sep=""))
} else {
  cor_mapa = ggraph(graph_corr, layout='nicely') +
    geom_edge_link(aes(edge_alpha = r, edge_width = r, color = r)) +
    scale_edge_colour_gradientn(limits = c(ccth, 1), colors = c(rgb(252,146,114,maxColorValue=255), "firebrick2")) +
    geom_node_point(fill = "white", size = 4,shape=21) + 
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_graph() +
    labs(title = paste("CROP - Correlation-based Removal of Multiplicities
Correlation network with correlation coefficient above ",ccth," and retention time window +- ",rtw, rtunit, sep=""))
}
  return(cor_mapa)
}

#######################################
# Function for multiplicities removal #
#######################################
mltpl_remove=function (data,det,name="",ccth,rtw,mcs,rt,maxrtw=NULL,funit,rtunit,mztab,mscsv,mztab_all){
  # Selection of clusters
  cor_tbl2 <- det
  if (mztab==TRUE){
    for (i in 1:nrow(cor_tbl2)) {
      cor_tbl2$x[i]=strsplit(cor_tbl2$x,"|",fixed=TRUE)[[i]][1]
      cor_tbl2$y[i]=strsplit(cor_tbl2$y,"|",fixed=TRUE)[[i]][1]
    }
  }
  
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
    
    mapa_pos=sort(which(colnames(data) %in% unique(c(cor_tbl2$x,cor_tbl2$y))))
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
      if (mztab==FALSE) {nms_mltpl = paste(nms_mltpl,colnames(data_mod)[clmn_del[k]], sep=" + ")
      } else nms_mltpl = paste(nms_mltpl,colnames(data_mod)[clmn_del[k]], sep="|")
    colnames(data_mod)[clmn_keep] = paste(colnames(data_mod)[clmn_keep],"*",sep="")
    data_mod[1,clmn_keep]=nms_mltpl
    clustn_short[i]=colnames(data_mod)[clmn_keep]
    clustn_long[i]=data_mod[1,clmn_keep]
    data_mod = data_mod[,-c(clmn_del)]
  }
 
  # Creating CROPPed tables
  if(mztab==TRUE) {
    data_mod2=mztab_save(data_mod,mscsv,name,mztab_all)
    write.table(data_mod2,paste(name,"_CROPped_ccth",ccth,"_rtw+-",rtw,".mzTab",sep=""), sep="\"", row.names=FALSE, col.names=FALSE, quote=FALSE) 
  } else {  
  
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
  }
  
  if (class(maxrtw) != "NULL") return(str_groups) #vector with stretched clusters coding
}

#####################################
# Function for loading mzTab format #
#####################################
mztab_load=function(mscsv=""){
  mztab_input= read.table(mscsv, sep="\"",blank.lines.skip = FALSE) 
  g_sml=grep("SML",mztab_input[,1])
  n=length(strsplit(as.character(mztab_input[g_sml[1],]),"\t",fixed=TRUE)[[1]])
  SML=matrix(ncol=n,nrow=length(g_sml))
  for (i in 1:length(g_sml))
    SML[i,]=strsplit(as.character(mztab_input[g_sml[i],]),"\t",fixed=TRUE)[[1]]
  
  g_smf=grep("SMF",mztab_input[,1])[-1]
  n=length(strsplit(as.character(mztab_input[g_smf[1],]),"\t",fixed=TRUE)[[1]])
  SMF=matrix(ncol=n,nrow=length(g_smf))
  for (i in 1:length(g_smf))
    SMF[i,]=strsplit(as.character(mztab_input[g_smf[i],]),"\t",fixed=TRUE)[[1]]
  
  data= SML[-1,grep("abundance_assay",SML[1,])] 
  colnames(data)= SML[1,grep("abundance_assay",SML[1,])]
  rownames(data)= SML[-1,grep("SML_ID",SML[1,])]
  ID=SML[-1,grep("SMF_ID",SML[1,])] #link to SMF table
  cl_names=SML[-1,grep("theoretical_neutral_mass",SML[1,])]
  
  SMF_pos=matrix(nrow=nrow(data),ncol=max(lengths(strsplit(ID,"|",fixed=TRUE))))
  for (i in 1:nrow(SMF_pos))
    for (j in 1:ncol(SMF_pos))
      SMF_pos[i,j] = strsplit(ID,"|",fixed=TRUE)[[i]][j]
  if (anyNA(SMF_pos[,1])==TRUE) stop("missing values in SMF_ID_REFS \n     cannot extract retention times from SMF table")
  SMF_pos2 = (as.numeric(unlist(SMF_pos)))
  SMF_pos= matrix(SMF_pos2, nrow=nrow(SMF_pos), ncol=ncol(SMF_pos))
  
  
  SMF_pos2=SMF_pos
  SMF_rt=suppressWarnings(as.numeric(SMF[-1,grep("retention_time",SMF[1,])[1]]))
  if (anyNA(SMF_rt)==TRUE) stop("missing values in retention_time_in_seconds/minutes column of SMF table \n     cannot compute correlations in RT windows")
  for (i in 1:nrow(SMF_pos2))
    for (j in 1:ncol(SMF_pos2))
      SMF_pos2[i,j]=SMF_rt[SMF_pos2[i,j]]
  
  SML_rt=rowMeans(SMF_pos2, na.rm=TRUE)
  data2 = suppressWarnings(matrix(as.numeric(unlist(data)), nrow=nrow(data), ncol=ncol(data)))
  rownames(data2)=rownames(data)
  colnames(data2)=colnames(data)
  if (anyNA(data2)==TRUE) {
    warning("missing values present in SML table \n     missing values treated as zeros when computing correlations")
    data2[is.na(data2)]=0               
  }
  data=cbind(rownames(data),SML_rt,data2)
  colnames(data)[1]="SML_ID"
  
  return(list(data=data,mztab_input=mztab_input,SML=SML,SMF_pos=SMF_pos,g_sml=g_sml,cl_names=cl_names))
}

####################################
# Function for saving mzTab format #
####################################
mztab_save=function(data,mscsv,name, mztab_all){
  data=t(data[-2,]) #data_mod from mltpl_remove function
  ID_SML=data[,1]
  data=data[,-1]
  
  SML_pos=matrix(nrow=length(ID_SML),ncol=max(lengths(strsplit(ID_SML,"|",fixed=TRUE))))
  for (i in 1:nrow(SML_pos))
    for (j in 1:ncol(SML_pos))
      SML_pos[i,j] = strsplit(ID_SML,"|",fixed=TRUE)[[i]][j]
  SML_pos= matrix(as.numeric(unlist(SML_pos)), nrow=nrow(SML_pos), ncol=ncol(SML_pos))
  
  new_pos=matrix(nrow=nrow(SML_pos),ncol=ncol(SML_pos)*ncol(mztab_all$SMF_pos))
  pos_con=numeric()
  for (i in 1:nrow(new_pos))
  {
    for (j in 1:ncol(SML_pos)) pos_con=c(pos_con,mztab_all$SMF_pos[SML_pos[i,j],])
    
    new_pos[i,]=pos_con
    pos_con=numeric()
  }  
  
  ID_SML_new=character()
  for(i in 1:nrow(new_pos))
    ID_SML_new[i]=paste0(na.omit(sort(unique(new_pos[i,]))),collapse="|")
  
  SML_new=mztab_all$SML[c(1,1+SML_pos[,1]),] #only rows with kept compounds
  SML_new[-1,grep("SMF_ID",mztab_all$SML[1,])]=ID_SML_new

  for (i in 1:nrow(SML_new)) SML_new[i,1]=stri_join(SML_new[i,],collapse="\t")
  
  mztab_output=as.character(mztab_all$mztab_input[,1])
  mztab_output[mztab_all$g_sml[1]:(mztab_all$g_sml[1]+nrow(SML_new)-1)]=SML_new[,1]
  mztab_output=mztab_output[-((mztab_all$g_sml[1]+nrow(SML_new)):(mztab_all$g_sml[1]+nrow(mztab_all$SML)-1))]
  mztab_output=as.data.frame(as.factor(mztab_output))
  
  return(mztab_output)
}

#########################################
# CROP - Function for the whole process #
#########################################
CROP=function(mscsv, name="", ccth = 0.75, rtw = 0.02, mcs=100, maxrtw = NULL, rtunit="min", funit="MW") {
  if (rtunit!="min" & rtunit!="s") stop("incorrect setting of rtunit   ['min' or 's']")
  mztab=FALSE
  
  if (length(grep(".csv",mscsv))>0) {
  data = read.csv2(mscsv,header=TRUE,dec=".")
  
  if (dim(data)[2]==1) data = read.csv2(mscsv,header=TRUE,dec=".",sep=",")
    
  if (dim(data)[2]==1) {
    stop ("only one column detected in the data \n   1: are you using a correct delimiter?   [';' or ','] \n   2: do you have first two columns with feature names \n      and retention time in the data?")
  } else if(dim(data)[2]<4) {
    stop("cannot compute correlations on 1 observation!")
  } else if(dim(data)[2]<7) {
    warning(paste("computing correlations only on ",dim(data)[2]-2," observations \n     consider using bigger sample size",sep=""))
  }
  
  if (dim(data)[1]==1) {
    stop("cannot compute correlations on 1 variable!")
  } else if ((dim(data)[1]<101) & (dim(data)[1]<dim(data)[2])) {
    warning(paste("computing correlations only on ",dim(data)[1]," variables \n     consider using more features \n  are you sure you have features in rows and samples in columns?",sep=""))
  } else if (dim(data)[1]<101) {
    warning(paste("computing correlations only on ",dim(data)[1]," variables \n     consider using more features",sep=""))
  } else if (dim(data)[1]<dim(data)[2]) {
    warning("are you sure you have features in rows and samples in columns?")
  }
  } else if (length(grep(".mzTab",mscsv))>0) {
    mztab=TRUE
    mztab_all=mztab_load(mscsv)
    data = mztab_all$data 
    }
  
  if (class(maxrtw) == "NULL") warning("you may have stretched clusters phenomenon in your data \n     which could not be highlighted in the correlation network \n  consider setting maxrtw for this purpose")
  
  rownames(data)=data[,1]
  data=data[,-1] 
  
  if(mztab==TRUE){
    data2 = data
    data = matrix(as.numeric(unlist(data)), nrow=nrow(data), ncol=ncol(data))
    colnames(data)=colnames(data2)
    rownames(data)=rownames(data2)
  }
    
  rt_vec=data[,1]  
  data=data[,-1] 
  data=t(data)
  
  if (mztab==FALSE) {rownames(data)=substring(rownames(data), 2)}
  
  names(rt_vec)=colnames(data)
  rt_vec=as.data.frame(t(rt_vec))
  data=as.data.frame(data)
  
  if (min(data)<0) warning("are you using non normalized data? \n     negative values detected")
  
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
  if(mztab==FALSE) {
    data_pro_mapu=mltpl_detect(data=data, ccth=ccth, rtw=rtw, rt=rt_vec)
  } else { 
    data_cl_n=data
    rt_vec_cl_n=rt_vec
    colnames(data_cl_n)=colnames(rt_vec_cl_n)=paste(colnames(data),mztab_all$cl_names,sep="|")
    data_pro_mapu=mltpl_detect(data=data_cl_n, ccth=ccth, rtw=rtw, rt=rt_vec_cl_n)
  }
  
  if (class(maxrtw)!="NULL") {
    sc_identif=mltpl_remove(data=data, det=data_pro_mapu, name=name, ccth=ccth, rtw=rtw, mcs=mcs, rt=rt_vec, maxrtw=maxrtw,rtunit=rtunit,funit=funit, mztab=mztab, mscsv=mscsv, mztab_all=mztab_all)
    cairo_pdf(paste(name=name,"_CROPped_ccth",ccth,"_rtw+-",rtw,"_correlation_network.pdf",sep=""),w=25,h=25)
    print(mltpl_depict(det=data_pro_mapu, ccth=ccth, rtw=rtw, stretched=sc_identif, rtunit=rtunit,mztab=mztab))
    dev.off()
  } else {
    mltpl_remove(data=data, det=data_pro_mapu, name=name, ccth=ccth, rtw=rtw, mcs=mcs, rt=rt_vec, rtunit=rtunit,funit=funit, mztab=mztab, mscsv=mscsv, mztab_all=mztab_all)
    cairo_pdf(paste(name=name,"_CROPped_ccth",ccth,"_rtw+-",rtw,"_correlation_network.pdf",sep=""),w=25,h=25)
    print(mltpl_depict(det=data_pro_mapu, ccth=ccth, rtw=rtw, rtunit=rtunit,mztab=mztab))
    dev.off()
  }
}

###########################################################################################################################################
###########################################################################################################################################
##################
# CROP IN ACTION #
##################
CROP(mscsv="example_data.csv", name="project1", ccth=0.75, rtw=0.02, maxrtw=2*0.02, rtunit="min")

###########################################################################################################################################
###########################################################################################################################################
#############################
# DESCRIPTION OF PARAMETERS #
#############################

###
# ccth ... threshold for correlation coefficient values 
#          default: ccth = 0.75
###
# funit ... units in which your features are assigned ("MW", "m/z" or anything else)
#           only for column headers in the csv outputs
#           ignore if using mzTab format
#           default: funit = "MW"
###
# maxrtw... maximal allowed RT window to color stretched clusters phenomenon
#           recommended to start with 2*rtw (higher number equals milder condition) 
#           default: maxrtw = NULL
###
# mcs ... maximal allowed cluster size
#         should not be set high; if default is not enough, rather consider setting ccth bigger and/or rtw smaller than changing mcs
#         default: mcs = 100
###
# mscsv ... filepath to an input mzTab or csv table of your MS data 
#           if using csv, make sure you have:
#              names of features in the first column
#              retention time of features in minutes in second column
#              samples in the rest of columns (including QCs)
#              column header as names of samples
#           if using mzTab, make sure you have:
#              values in SML table after imputation of zeros
#              no missing values in columns with SML_ID and SMF_ID_REFS of SML table
#              no missing values in columns with SMF_ID and retention_time_in_seconds/minutes of SMF table
###
# name ... a note which will be in names of all files
###
# rtunit ... units of retention time you ar using in your data ("min" or "s")
#            default: rtunit = "min"
###
# rtw ... retention time window where +-rtw will be considered
#         default: rtw = 0.02



