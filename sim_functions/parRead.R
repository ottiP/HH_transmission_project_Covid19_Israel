parRead <- function(N){
  summary_tables <- readRDS(paste('../intermediate_data/summary_table_',N,'iter_SIM.rds'))
  output <- model.run(summary_tables)
  return(output)
}