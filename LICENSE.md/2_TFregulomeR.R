library(TFregulomeR)
all_record <- TFBSBrowser()

# download all peaks with motif
for (i in all_record$ID){
peak_i = loadPeaks(id=i, includeMotifOnly=T)
write.table(peak_i, paste0(i,'_peaks_with_motif.txt'), sep='\t', col.names=F, row.names=F, quote=F)}

# Downlload all peaks
for (i in all_record$ID){
peak_i = loadPeaks(id=i, includeMotifOnly=F)
write.table(peak_i, paste0(i,'_peaks_all.txt'), sep='\t', col.names=F, row.names=F, quote=F)}

