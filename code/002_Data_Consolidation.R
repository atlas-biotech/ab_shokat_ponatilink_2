# 002 Data Consolidation #

# Load dependencies #
library(dplyr)

# Read in merged dataframes #
df.twist = read.csv('code/Refs/Twist_ABL_All_Data.csv', header=T)
seq.bl.a = read.csv('merged_data/Full_SH23K_K_BL_A/TileSeq_K562_full_kinase_SH2_3_baseline_D0_A.csv', header=T)
seq.a.1000.a = read.csv('merged_data/Full_SH23K_KRA_Asc_1000_A/TileSeq_K562_full_kinase_SH2_3_KRA_Asc_1000_A.csv', header=T)
seq.pl2.100.a = read.csv('merged_data/Full_SH23K_PL2_100_A/TileSeq_K562_full_kinase_SH2_3_PL2_100_A.csv', header=T)
seq.p.30.a = read.csv('merged_data/Full_SH23K_P30_D21_A/TileSeq_K562_full_kinase_SH2_3_P30_D21_A.csv', header=T)
seq.bl.b = read.csv('merged_data/Full_SH23K_K_BL_B/TileSeq_K562_full_kinase_SH2_3_KE1_BL_B.csv', header=T)
seq.a.1000.b = read.csv('merged_data/Full_SH23K_KRA_Asc_1000_B/TileSeq_K562_full_kinase_SH2_3_KRA_Asc_1000_B.csv', header=T)
seq.pl2.100.b = read.csv('merged_data/Full_SH23K_PL2_100_B/TileSeq_K562_full_kinase_SH2_3_PL2_100_B.csv', header=T)
seq.p.30.b = read.csv('merged_data/Full_SH23K_P30_D21_B/TileSeq_K562_full_kinase_SH2_3_P30_D21_B.csv', header=T)

df.names = grep("^seq\\.*", ls(), value = TRUE)
df.list = mget(df.names)

df.list.proc = lapply(df.list, function(df) {
  df %>% mutate(mut = paste0(ref_aa,protein_start,alt_aa)) %>% 
                # freq = ct/depth) %>% 
    filter(ct > 10, alt_aa != 'O') %>% 
    select(mut,alt_codon,ct,depth) %>% 
    group_by(mut, alt_codon) %>%
    summarise(
      ct = sum(ct),                # Sum ct values
      depth = mean(depth),         # Average depth values
      .groups = "drop") %>% 
    mutate(freq = ct/depth)
})


df.out = df.twist %>% 
  select(-AA.Position) %>%
  dplyr::rename(alt_codon = variant_codon, alt_aa = variant_aa)

for(i in 1:length(df.list.proc)) {
  df = df.list.proc[[i]]
  input = df.names[i]
  condition =  sub("^seq", "", input)
  newnames = c('ct','depth','freq')
  colnames(df)[colnames(df) %in% newnames] = paste0(colnames(df)[colnames(df) %in% newnames], condition)
  df.out = left_join(df.out,df,by = c('mut','alt_codon'))
}

# pseudocount.a = min(df.out$freq.bl.a, na.rm=T)
# pseudocount.b = min(df.out$freq.bl.b, na.rm=T)

# df.out.pseudo = df.out %>%
#   mutate(
#     freq.bl.a = if_else(
#       is.na(freq.bl.a) & (!is.na(freq.a.1000.a) | !is.na(freq.pl2.100.a) | !is.na(freq.p.30.a)),
#       pseudocount.a,
#       freq.bl.a
#     ),
#     freq.bl.b = if_else(
#       is.na(freq.bl.b) & (!is.na(freq.a.1000.b) | !is.na(freq.pl2.100.b) | !is.na(freq.p.30.b)),
#       pseudocount.b,
#       freq.bl.b
#     ))

df.out.proc = df.out %>% 
  rowwise() %>% 
  mutate(
    freq.bl = mean(c(freq.bl.a,freq.bl.b),na.rm=T),
    freq.a.1000 = mean(c(freq.a.1000.a,freq.a.1000.b),na.rm=F),
    freq.pl2.100 = mean(c(freq.pl2.100.a,freq.pl2.100.b),na.rm=F),
    freq.p.30 = mean(c(freq.p.30.a,freq.p.30.b),na.rm=F),
  ) %>% 
  ungroup()

df.out.l2fc = df.out.proc %>% 
  mutate(
    l2fc.a.1000 = log2(freq.a.1000/freq.bl),
    l2fc.pl2.100 = log2(freq.pl2.100/freq.bl),
    l2fc.p.30 = log2(freq.p.30/freq.bl)
  )

write.csv(x = df.out.l2fc, file = 'processed_data/pl2_master_datafile.cvs', row.names = F)

df.shokat = df.out.proc %>%
  select(
    mut,
    wt_codon,
    wt_aa,
    alt_codon,
    alt_aa,
    variant_proportion,
    pos.aa,
    snp,
    ct.bl.a,
    depth.bl.a,
    freq.bl.a,
    ct.bl.b,
    depth.bl.b,
    freq.bl.b,
  )

write.csv(x = df.shokat, file = 'processed_data/ABL_baseline_master_file.csv', row.names = F)
