library(tidyverse)
library(data.table)

#### Sep 06, 2023 ####

# counts = matrix of hom-ref, het, hom-alt with snp name as rowname
# calculates -log10p
hwe_delta <- function(counts){
  
  outmat = as.data.table(bind_cols(ID=rownames(counts),n=rowSums(counts),counts))
  names(outmat) = c("ID","n","aa_ct","ab_ct","bb_ct")
  
  outmat[,aa_fq:=aa_ct/n]
  outmat[,ab_fq:=ab_ct/n]
  outmat[,bb_fq:=bb_ct/n]
  outmat[,a_fq:=aa_fq+ab_fq/2]
  outmat[,b_fq:=bb_fq+ab_fq/2]
  outmat[,delta:=aa_fq-a_fq^2]
  outmat[,se:=sqrt((1/n) * a_fq^2 * (1-a_fq)^2)]
  outmat[,stat:=(delta/se)^2]
  outmat[,log10p:=-pchisq(stat, df =1, lower.tail = F, log.p = T)/log(10)]
  
  return(outmat)
}

# dese = named list of hwd matrices (delta, se, n, stat) with groups as colnames and snps as rownames
metas <- function(dese){
  
  invar_param = (1/dese$se)^2
  
  # inverse variance
  
  outmet = data.table(ID = dese$ID, invar_se = sqrt(1/rowSums(invar_param)), invar_beta = rowSums(dese$delta*invar_param)/rowSums(invar_param))
  outmet[,invar_stat:=invar_beta/invar_se]
  outmet[,invar_log10p:=-(log(2)+pnorm(abs(invar_stat),lower.tail = F, log.p=T))/log(10)]
  
  # sample size
  
  outmet[,samsi_stat:=rowSums((dese$delta/dese$se) *sqrt(dese$n))/ sqrt(rowSums(dese$n))]
  outmet[,samsi_log10p:=-(log(2)+pnorm(abs(samsi_stat), lower.tail = F, log.p=T))/log(10)]
  
  # omnibus 
  
  outmet[,omni_stat:=rowSums(dese$stat)]
  outmet[,omni_log10p:=-pchisq(omni_stat, df=ncol(dese$se), lower.tail = F, log.p = T)/log(10)]
  
  # heterogeneity
  
  outmet[,hegen_stat:=omni_stat-1/rowSums(invar_param)*rowSums(dese$delta*invar_param)^2]
  outmet[,hegen_log10p:=-pchisq(hegen_stat, df=(ncol(dese$se)-1), lower.tail = F, log.p = T)/log(10)]
  
  
  return(outmet)
}

# order should be female=1, male=2
autohwd1df <- function(feme){
  outwld = data.table(ID = feme$ID)
  
  pF=feme$b_fq[[1]]
  pM=feme$b_fq[[2]]
  nF=feme$n[[1]]
  nM=feme$n[[2]]
  delta.F=feme$delta[[1]]
  delta.M=feme$delta[[2]]
  
  outwld[,sdmaf_stat:=(pM-pF)^2/(1/(2*nM)*(pM*(1-pM)+delta.M)+1/(2*nF)*(pF*(1-pF)+delta.F))]
  outwld[,sdmaf_log10p:=-pchisq(sdmaf_stat,df=1,lower.tail = F,log.p=T)/log(10)]
  
  maf.F=apply(matrix(c(a1f=feme$b_fq[[1]], a2f=1-feme$b_fq[[1]]), ncol=2),1,min)
  maf.M=apply(matrix(c(a1f=feme$b_fq[[2]], a2f=1-feme$b_fq[[2]]), ncol=2),1,min)
  
  outwld[,sdmaf_diff:=maf.F-maf.M] # order matters here
  
  return(outwld)
}

allmet <- function(hwe_com){
  
  deltas = map(hwe_com, dplyr::select, ID, delta, se, stat, n, b_fq) %>% map(tibble::column_to_rownames, "ID") %>% map(as.matrix)
  
  commsnp = purrr::reduce(sapply(map(deltas,~.[complete.cases(.),]),rownames),intersect)
  melta=NULL
  waldo=NULL
  if(length(commsnp)>0) {
    message("Performing meta-analysis on ", length(commsnp), " snps with non-missing p-value in all groups")
    
    deltasub = NULL
    deltasub[["ID"]] = commsnp
    for(i in colnames(deltas[[1]])){
      deltasub[[i]] = as.data.frame(map(deltas,~.[commsnp,i]))
    }
    melta = metas(deltasub)
    if(length(hwe_com)==2 & grepl("\\.male",names(hwe_com)[2]) & grepl("female",names(hwe_com)[1])){
      message("Adding sdMAF-analysis")
      waldo = autohwd1df(deltasub)
    }
  }
  return(list(melta=melta,waldo=waldo))
}


## combine ##
# groups = named list of matrices of hom-ref, het, hom-alt with snp name as rowname
hwe_meta_bygroup <- function(groups, clusters){
  message("Num snps per group : ")
  for(i in names(groups)) message(i, " = ", nrow(groups[[i]]))
  
  message("===============")
  message("Calculating HWD for each group separately")
  message("===============")
  hwe_out = map(groups, hwe_delta)
  met_out = NULL
  if(is.null(clusters)) met_out = allmet(hwe_out)
  for(g in names(clusters)){
    message("Group membership : ",g)
    met_out[[g]] = allmet(hwe_out[clusters[[g]]])
  }
  
  met_outlist = purrr::list_transpose(met_out)
  
  return(list(hwe = hwe_out, meta = rlist::list.clean(met_outlist$melta), sdmaf = rlist::list.clean(met_outlist$waldo)))
}
