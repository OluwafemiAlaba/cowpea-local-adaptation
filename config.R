# Run scripts from the repository root. Keep the extracted Dryad directory intact.
project <- path.expand(Sys.getenv('COWPEA_ROOT', '~/data/cowpealandrace_460'))
inputs <- file.path(project, 'manuscript_reanalysis')
outputs <- path.expand(Sys.getenv('COWPEA_OUTPUT', file.path(project, 'reproduced')))
stopifnot(dir.exists(inputs), startsWith(project, '/'), startsWith(outputs, '/'))
folders <- c(structure='genetic_structure_verified_460', rda='rda_verified_460',
             egwas='egwas_verified_460', go='go_window20kb', gf='GF_validation_design_review')
input <- function(stage, ...) file.path(inputs, folders[[stage]], ...)
output <- function(stage, ...) {
  p <- file.path(outputs, stage, ...)
  dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)
  p
}
# Explicit opt-in lets a single downstream step use the deposited results.
result <- function(stage, ...) {
  p <- output(stage, ...)
  if (!file.exists(p) && Sys.getenv('COWPEA_USE_DEPOSITED') == '1') p <- input(stage, ...)
  if (!file.exists(p)) stop('Required preceding result: ', p)
  p
}
write_tsv <- function(x, p) write.table(x, p, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
