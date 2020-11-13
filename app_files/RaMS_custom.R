# filename <- "G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS\\190715_Poo_TruePooFK180310_DDApos20.mzML"
# filename <- "G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS\\180205_Poo_TruePooPos_dda3.mzML"


grabMzmlMetadata <- function(xml_data){
  compr_xpath <- paste0('//d1:cvParam[@accession="MS:1000574"]|',
                        '//d1:cvParam[@accession="MS:1000576"]')
  compr_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, compr_xpath), "name")
  compr <- switch(compr_type,
                  `zlib compression`="gzip",
                  `no compression`="none",
                  `none`="none")
  
  mz_precision_xpath <- '//d1:cvParam[@accession="MS:1000523"]'
  mz_bit_type <- xml_attr(xml_find_first(xml_data, mz_precision_xpath), "name")
  mz_precision <- sub(mz_bit_type, pattern = "-bit float", replacement = "")
  mz_precision <- as.numeric(mz_precision)/8
  
  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, int_bit_xpath), "name")
  int_precision <- sub(int_bit_type, pattern = "-bit float", replacement = "")
  int_precision <- as.numeric(int_precision)/8
  
  ms2_ce_xpath <- '//d1:cvParam[@accession="MS:1000045"]'
  ms2_ce_val <- xml2::xml_attr(xml2::xml_find_first(xml_data, ms2_ce_xpath), "value")
  ms2_ce_val <- as.numeric(ms2_ce_val)
  
  list(compression=compr, mz_precision=mz_precision, int_precision=int_precision,
       ms2_ce=ms2_ce_val)
}

grabMzmlData <- function(filename){
  xml_data <- xml2::read_xml(filename)
  
  file_metadata <- grabMzmlMetadata(xml_data)
  
  ms1_nodes <- xml2::xml_find_all(
    xml_data, '//d1:cvParam[@name="ms level"][@value="1"]/parent::d1:spectrum'
  )
  
  rt_vals <- grabSpectraRt(ms1_nodes)
  mz_vals <- grabSpectraMz(ms1_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms1_nodes, file_metadata)
  
  ms1_data <- data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
                         mz=unlist(mz_vals), int=unlist(int_vals))
  
  ms2_xpath <- '//d1:cvParam[@name="ms level"][@value="2"]/parent::d1:spectrum'
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  
  rt_vals <- grabSpectraRt(ms2_nodes)
  premz_vals <- grabSpectraPremz(ms2_nodes)
  mz_vals <- grabSpectraMz(ms2_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms2_nodes, file_metadata)
  
  ms2_data <- data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
                         premz=rep(premz_vals, sapply(mz_vals, length)),
                         fragmz=unlist(mz_vals), int=unlist(int_vals),
                         voltage=file_metadata$ms2_ce)
  
  list(ms1_data, ms2_data)
}

pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}

grabSpectraRt <- function(xml_nodes){
  rt_xpath <- 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(xml_nodes, rt_xpath)
  as.numeric(xml2::xml_attr(rt_nodes, "value"))
}

grabSpectraPremz <- function(xml_nodes){
  premz_xpath <- paste0('d1:precursorList/d1:precursor/d1:selectedIonList',
                        '/d1:selectedIon/d1:cvParam[@name="selected ion m/z"]')
  premz_nodes <- xml2::xml_find_all(xml_nodes, premz_xpath)
  as.numeric(xml2::xml_attr(premz_nodes, "value"))
}

grabSpectraMz <- function(xml_nodes, file_metadata){
  mz_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[1]/d1:binary'
  mz_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, mz_xpath))
  lapply(mz_vals, function(binary){
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$mz_precision,
                            size = file_metadata$mz_precision)
  })
}

grabSpectraInt <- function(xml_nodes, file_metadata){
  int_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[2]/d1:binary'
  int_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, int_xpath))
  int_vals <- lapply(int_vals, function(binary){
    decoded_binary <- base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$int_precision,
                            size = file_metadata$int_precision)
  })
}

