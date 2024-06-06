# Read the two files that we want to merge
occ1 <- read.csv(paste0("../../data/ptealc/ebd_ptealc_breeding_spain_zf_part1_variables.csv"))
occ2 <- read.csv(paste0("../../data/ptealc/ebd_ptealc_breeding_spain_zf_part2_variables.csv"))

# Merge both files
occ_bind <- rbind(occ1, occ2)

# Delete the last column named ".geo", it is a GeoJSON point we won't need it
# Also delete system.index column
columns_to_remove <- c(".geo", "system.index")
occ_bind <- occ_bind[, !names(occ_bind) %in% columns_to_remove]

# Save the merged file
write.csv(occ_bind, file = paste0("../../data/ptealc/ebd_ptealc_breeding_spain_variables.csv"), row.names = FALSE, quote=FALSE)