time /home/hr/openmsall/builds/openms_trunk2/bin/OpenSwathWorkflow "-in" /tmp/qq*.mzML  "-tr" /media/doc/project_data/../tmp/openms_issues/Hannes_slow_osw/sublib2X.tsv  -tr_irt   "-out_chrom" "/tmp/OpenSwathWorkflow_1.chrom.mzML.tmp" -out_tsv /tmp/out.tsv  -readOptions workingInMemory  -force   -Scoring:TransitionGroupPicker:compute_peak_quality true -Scoring:TransitionGroupPicker:minimal_quality -9999  -use_ms1_traces -rt_norm /tmp/test.trafoXML -min_upper_edge_dist 1  -Scoring:stop_report_after_feature -1 
source ~/pyenv/py27_bleeding/bin/activate
cp /tmp/out.tsv /tmp/out.all3.tsv
pyprophet /tmp/out.all3.tsv  --ignore.invalid_score_columns

# Jul 10
# then do this in R:
      /*

df = read.csv("/tmp/out_with_dscore_filtered.csv", sep="\t")
df = read.csv("/tmp/out.all2_with_dscore_filtered.csv", sep="\t")
df = read.csv("/tmp/out.all3_with_dscore_filtered.csv", sep="\t")
df = read.csv("/media/doc/tmp/openms_issues/output_xall_scored.tsv", sep="\t")
df = fread("/media/doc/tmp/openms_issues/output_xall_scored.tsv", sep="\t")
names(df)
df = subset(df, d_score > -6.0)

df_filter = subset(df, q_value < 0.3)
df_orig 

# check if we really need the top5 peaks all the time ...
plot(df$d_score, df$main_var_xx_swath_prelim_score, cex=0.1); abline(h=0.5); abline(v=1.9)
plot(df$d_score, df$peak_group_rank, cex=0.1); abline(h=0.5); abline(v=1.9) 

# plot(df$d_score, df$peak_group_rank, cex=0.1); abline(h=0.5); abline(v=1.9) 

hist(df$main_var_xx_swath_prelim_score, cex=0.1)
quantile(df$main_var_xx_swath_prelim_score)
#          0%         25%         50%         75%        100%
# -6.33113000 -0.08331215  0.39149400  0.93274925  7.68394000

plot(df$d_score, df$main_var_xx_swath_prelim_score, cex=0.1); abline(h=0.5); abline(v=1.9)
plot(df$d_score, df$initialPeakQuality, cex=0.1); abline(v=2.9); abline(h=-0.5)
plot(df$d_score, df$init_shape_score, cex=0.1)
plot(df$d_score, df$init_coel_score, cex=0.1)
plot(df$d_score, df$var_xcorr_coelution, cex=0.1); abline(v=2.9)
plot(df$d_score, df$var_xcorr_shape, cex=0.1); abline(v=2.9)
plot(df$d_score, df$var_library_sangle, cex=0.1); abline(v=2.9)
plot(df$d_score, df$var_library_rmsd, cex=0.1); abline(v=2.9)
plot(df$d_score, df$init_missing_peaks, cex=0.1); abline(v=2.9)
dim(subset(df, m_score < 0.01)) # 66
dim(subset(df, m_score < 0.1)) # 138

df10 = (subset(df, m_score < 0.1)) # 138

dim(subset(df, main_var_xx_swath_prelim_score < 0.5 & d_score > 1))
dim(subset(df, main_var_xx_swath_prelim_score < 0.5 & d_score > 2)) # 140
dim(subset(df, main_var_xx_swath_prelim_score < 0.4 & d_score > 1))
dim(subset(df, main_var_xx_swath_prelim_score < 0.4 & d_score > 1.5)) # 511
dim(subset(df, main_var_xx_swath_prelim_score < 0.4 & d_score > 2)) # 100

min(subset(df, d_score < 2)$q_value)   # 0.2
min(subset(df, d_score < 1.5)$q_value) # 0.3
min(subset(df, d_score < 1)$q_value)   # 0.44


hist(df$delta_rt)
hist(subset(df, d_score > 0.05)$delta_rt)
hist(subset(df, d_score > 0.05)$delta_rt)
plot(df10$RT, df10$delta_rt)
plot(df10$RT, df10$norm_RT)
plot(df10$RT, df10$assay_rt)

dim(subset(df, initialPeakQuality < -0.5 & d_score > 2.9))

plot(df$init_shape_score, df$var_xcorr_coelution, cex=0.1); abline(v=2.9)

min((subset(df, q_value < 0.1)$d_score) ) # 2.9

dim(subset(df, main_var_xx_swath_prelim_score < 0)) # 50k
dim(subset(df, main_var_xx_swath_prelim_score > 0)) # 12k

df2 = subset(df, init_coel_score > 0.2 & init_shape_score  < 0.4)
plot(df2$d_score, df2$init_missing_peaks, cex=0.1); abline(v=2.9)

library(MASS)
df_n = df
df_n$pre = df_n$d_score > 2.9
l = lda(pre ~ init_coel_score + init_shape_score + init_missing_peaks, df_n)
k = predict(l, df_n)
df_n$lda = predict(l, df_n)$x
plot(df_n$d_score, df_n$lda, cex=0.1); abline(v=2.9)
dim(subset(df_n, lda < 0.5 & d_score > 2.9)) # 6 FP
dim(subset(df_n, lda < 0.5 )) / dim(df_n) # removes 70% of the data

l = lda(pre ~ var_xcorr_coelution + var_xcorr_shape, df_n)
k = predict(l, df_n)
df_n$lda = predict(l, df_n)$x
plot(df_n$d_score, df_n$lda, cex=0.1); abline(v=2.6)

#                             LD1
# var_xcorr_coelution -0.1989225
# var_xcorr_shape      2.7652960


#     LD1
#     var_xcorr_coelution -0.1594062
#     var_xcorr_shape      1.5442088
#     var_library_sangle  -2.1177013


dim(subset(df_n, lda < 0.5 & d_score > 2.6)) # 2000 FN
dim(subset(df_n, lda < 0.5 )) / dim(df_n) # removes 73% of the data

dim(subset(df_n, lda < 0.0 & d_score > 2.6)) # 800 FN
dim(subset(df_n, lda < 0.0 )) / dim(df_n) # removes 73% of the data

dim(subset(df_n, lda < -0.5 & d_score > 2.6)) # 300 FP
dim(subset(df_n, lda < -0.5 )) / dim(df_n) # removes 35% of the data

dim(subset(df_n, lda < -1.0 & d_score > 2.6)) # 37 FP
dim(subset(df_n, lda < -1.0 )) / dim(df_n) # removes 15% of the data

dim(subset(df_n, lda < -1.5 & d_score > 2.6)) # 2 FP
dim(subset(df_n, lda < -1.5 )) / dim(df_n) # removes 3% of the data

l = lda(pre ~ var_xcorr_coelution + var_xcorr_shape + var_library_sangle, df_n)
k = predict(l, df_n)
df_n$lda = predict(l, df_n)$x
plot(df_n$d_score, df_n$lda, cex=0.1); abline(v=2.6)


#     LD1
#     var_xcorr_coelution -0.1594062
#     var_xcorr_shape      1.5442088
#     var_library_sangle  -2.1177013


dim(subset(df_n, lda < 0.5 & d_score > 2.6)) # 231 FN
dim(subset(df_n, lda < 0.5 )) / dim(df_n) # removes 69% of the data

dim(subset(df_n, lda < 0.0 & d_score > 2.6)) # 80 FN
dim(subset(df_n, lda < 0.0 )) / dim(df_n) # removes 50% of the data

dim(subset(df_n, lda < -0.5 & d_score > 2.6)) # 40 FP
dim(subset(df_n, lda < -0.5 )) / dim(df_n) # removes 35% of the data

dim(subset(df_n, lda < -1.0 & d_score > 2.6)) # 19 FP
dim(subset(df_n, lda < -1.0 )) / dim(df_n) # removes 15% of the data

dim(subset(df_n, lda < -1.5 & d_score > 2.6)) # 10 FP
dim(subset(df_n, lda < -1.5 )) / dim(df_n) # removes 5% of the data



*/
