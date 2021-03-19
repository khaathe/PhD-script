biotools.ens.to.symbol <- function(ensembl_ids){
  url = "https://biotools.fr/human/ensembl_symbol_converter/"
  ids_json <- toJSON(ensembl_ids)
  body <- list(api=1, ids=ids_json)
  r <- POST(url, body = body)
  output <- fromJSON( httr::content(r, "text"), flatten=TRUE)
  output <- output[-which(base::sapply(output, FUN = is.null))] # removing null values
  df <- as.data.frame(output)
  df <- data.frame(ensembl_ids = colnames(genes_name), symbol = t(output[1,]))
  return(df)
}

annotdbi.ens.to.symbol <- function(ens){
  df <- select(x = org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", columns = c("SYMBOL"))
  return(df)
}