# The following code have been adopted from https://github.com/Pinlyu3/IReNA-v2 
# and changed to suit our data

### Read fragments ####
Read_fragment_to_GR <- function(x){
  GR <- GRanges(seqnames=x$V1,ranges=IRanges(start=x$V2,end=x$V3),index=x$V4)
  return(GR)
}


Read_peak_to_GR <- function(x){
  GR <- GRanges(seqnames=x$V1,ranges=IRanges(start=x$V2,end=x$V3))
  return(GR)
}


#####

Check_normalized_Signal <- function(file,savefile){
  library(rtracklayer)
  temp_bw = import.bw(file)
  print(summary(width(temp_bw)))
  print(summary(temp_bw$score))
  setwd('/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/All_normalized')
  saveRDS(temp_bw,file=savefile)
}


#####
Calculate_signal_bw <- function(GR,Signal){
  ####
  res = findOverlaps(GR,Signal)
  res_score = Signal$score[subjectHits(res)]
  ####
  res_score_merge = tapply(res_score,queryHits(res),sum)
  ####
  GR$score = 0
  ####
  GR$score[as.numeric(names(res_score_merge))] = as.numeric(res_score_merge)
  ####
  return(GR$score)
}


Calculate_footprint <- function(footprint_GR,Signal){
  width = width(footprint_GR)
  #####
  left_s = start(footprint_GR)-width*3
  left_e = start(footprint_GR)-1
  #####
  left_GR = footprint_GR
  start(left_GR) = left_s
  end(left_GR) = left_e
  ######
  right_s = end(footprint_GR)+1
  right_e = end(footprint_GR)+width*3
  ######
  right_GR = footprint_GR
  start(right_GR) = right_s
  end(right_GR) = right_e
  ######
  center_GR = footprint_GR
  ######
  ######
  print('Calculate_footprint_center')
  footprint_GR$center = Calculate_signal_bw(center_GR,Signal)
  print('Calculate_footprint_left')
  footprint_GR$left = Calculate_signal_bw(left_GR,Signal)
  print('Calculate_footprint_right')
  footprint_GR$right = Calculate_signal_bw(right_GR,Signal)
  ######
  return(footprint_GR)
}


Calculate_footprint_celltypes <- function(footprint_GR,Signal,TFs,out_all_ext){
  motifs_need = out_all_ext$Motif[which(out_all_ext$TFs %in% TFs == T)]
  ####
  print(length(motifs_need))
  ####
  footprint_GR_cl = footprint_GR[which(footprint_GR$motifs %in% motifs_need == T)]
  ####
  footprint_GR_cl_out = Calculate_footprint(footprint_GR_cl,Signal)
  ####
  return(footprint_GR_cl_out)
}


Filter_footprints <- function(footprint_GR_cl_out,delta=0.1){
  ####
  left_delta = (footprint_GR_cl_out$left)/3 - footprint_GR_cl_out$center
  right_delta = (footprint_GR_cl_out$right)/3 - footprint_GR_cl_out$center
  ####
  k = which(footprint_GR_cl_out$left > delta*3 & footprint_GR_cl_out$right > delta*3 & footprint_GR_cl_out$center < -delta)
  ####
  footprint_GR_cl_out_filtered = footprint_GR_cl_out[k]
  ####
  return(footprint_GR_cl_out_filtered)
}