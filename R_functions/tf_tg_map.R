# Find overlaps of promoter peaks and TF-motif binding sites
tf_tg_map <- function(footprint_res_cl, Out_Dat_GR){
 
  #### merge Out_Dat with footprint_res_cl
  #### remove their not overlapped regions 
  k1 = which(countOverlaps(Out_Dat_GR, footprint_res_cl) >0)
  k2 = which(countOverlaps(footprint_res_cl, Out_Dat_GR) >0)
  
  #### keep overlapped
  Out_Dat_GR_Overlap = Out_Dat_GR[k1]
  footprint_res_cl_Overlap = footprint_res_cl[k2]
  
  #### merge Out_Dat_GR_Overlap and footprint_res_cl_Overlap
  Res_find = data.frame(findOverlaps(footprint_res_cl_Overlap, Out_Dat_GR_Overlap))
  
  ####
  Res_find$Motifs = footprint_res_cl_Overlap$motifs[Res_find$queryHits]
  Res_find$footprint = as.character(footprint_res_cl_Overlap)[Res_find$queryHits]
  
  ####
  Res_find$Target = Out_Dat_GR_Overlap$Target[Res_find$subjectHits]
  Res_find$peaks = as.character(Out_Dat_GR_Overlap)[Res_find$subjectHits]
  
  return(Res_find)
}



