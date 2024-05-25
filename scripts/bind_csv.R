# Define the list of the studied species and their abbreviated names
especies <- list(
  #c("Chersophilus duponti", "chedup"),   # Alondra ricotí / Dupont’s Lark
  #c("Circus pygargus", "cirpyg"),        # Aguilucho cenizo / Montagu's Harrier
  #c("Falco naumanni", "falnau"),         # Cernícalo primilla / Lesser Kestrel
  #c("Otis tarda", "ottard"),             # Avutarda euroasiática / Great Bustard
  #c("Pterocles alchata", "ptealc"),      # Ganga ibérica / Pin-tailed Sandgrouse
  #c("Pterocles orientalis", "pteori"),   # Ganga ortega / Black-bellied Sandgrouse
  c("Tetrax tetrax", "tettet")            # Sisón común / Little Bustard
)

# Loop over each species files to merge their two parts
for (i in seq_along(especies)) {

  # Read the two files that we want to merge
  df1 <- read.csv(paste0("data/",especies[[i]][2],"/ebd_",especies[[i]][2],"_breeding_spain_zf_part1_variables.csv"))
  df2 <- read.csv(paste0("data/",especies[[i]][2],"/ebd_",especies[[i]][2],"_breeding_spain_zf_part2_variables.csv"))
  
  # Merge both files
  df_bind <- rbind(df1, df2)
  
  # Delete the last column named ".geo", it is a GeoJSON point we won't need it
  df_bind <- df_bind[, -ncol(df_bind)]
  
  # Save the merged file
  write.csv(df_bind, file = paste0("data/",especies[[i]][2],"/ebd_", especies[[i]][2], "_breeding_spain_variables.csv"), row.names = FALSE)
}