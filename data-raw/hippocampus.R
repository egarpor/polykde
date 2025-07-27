
load("spokes.rda")

hippocampus <- list(
  "base" = base,
  "bdry" = bdry,
  "dirs" = dirs,
  "rads" = rads,
  "ids" = ids,
  "ids_labs" = ids_labs
)

# Save tables
save(list = "hippocampus",
     file = "../data/hippocampus.rda", compress = "xz")
