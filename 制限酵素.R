setwd("~/20250610_newtarget")

library(Biostrings)

# 1. IUPAC展開関数（GTMKAC → GTAAAC, GTACTCなど）
expand_iupac <- function(pattern) {
  chars <- strsplit(pattern, "")[[1]]
  options <- lapply(chars, function(ch) {
    bases <- Biostrings::IUPAC_CODE_MAP[[ch]]
    strsplit(bases, "")[[1]]
  })
  expand.grid(options, stringsAsFactors = FALSE) |>
    apply(1, paste0, collapse = "")
}

# 2. 制限酵素認識配列
restriction_enzymes <- list(
  AccI    = "GTMKAC",
  AflII   = "CTTAAG",
  AluI    = "AGCT",
  ApaI    = "GGGCCC",
  BamHI   = "GGATCC",
  BglII   = "AGATCT",
  BlnI    = "CCTAGG",
  Bsp1407I = "TGTACA",
  Bst1107I = "GTATAC",
  ClaI    = "ATCGAT",
  DpnI    = "GATC",
  DraI    = "TTTAAA",
  EcoO1091I = "RGGNCCY",
  EcoRI   = "GAATTC",
  EcoRV   = "GATATC",
  EcoT14I = "CCWGG",
  EcoT22I = "ATGCAT",
  HindIII = "AAGCTT",
  HpaI    = "GTTAAC",
  KpnI    = "GGTACC",
  MboI    = "GATC",
  MluI    = "ACGCGT",
  NaeI    = "GCCGGC",
  NcoI    = "CCATGG",
  NdeI    = "CATATG",
  NheI    = "GCTAGC",
  NotI    = "GCGGCCGC",
  PstI    = "CTGCAG",
  SacI    = "GAGCTC",
  SacII   = "CCGCGG",
  SalI    = "GTCGAC",
  Sau3AI  = "GATC",
  ScaI    = "AGTACT",
  SfiI    = "GGCCNNNNNGGCC",
  SmaI    = "CCCGGG",
  SpeI    = "ACTAGT",
  SphI    = "GCATGC",
  SspI    = "AATATT",
  StuI    = "AGGCCT",
  XbaI    = "TCTAGA",
  XhoI    = "CTCGAG",
  AscI    = "GGCGCGCC",
  AvrII   = "CCTAGG",
  AhdI    = "GACNNNNNGTC",
  AsiSI   = "GCGATCGC",
  AgeI    = "ACCGGT",
  BbsI    = "GAAGAC",
  BsaI    = "GGTCTC",
  BstI    = "TCGA",
  BstBI   = "TTCGAA",
  BsmBI   = "CGTCTC",
  BsrGI   = "TGTACA",
  BbeI    = "GGCGCC",
  BspQI   = "GCTCTTC",
  BslII   = "CCNNNNNNNGG",
  BcoDI   = "GTCTC",
  BsiWI   = "CGTACG",
  DraIII  = "CACNNNGTG",
  EcoNI   = "CCTNNNNNAGG",
  NspI    = "RCATGY",
  MseI    = "TTAA",
  MfeI    = "CAATTG",
  PciI    = "ACATGT",
  PaqCI   = "CCTT",
  SbfI    = "CCTGCAGG",
  SfoI    = "GGCGCC",
  SrfI    = "GCCCGGGC",
  SapI    = "GCTCTTC",
  XmnI    = "GAANNNNTTC",
  BanI    = "GGYRCC"
)

# 3. FASTA読み込み＋結合
fasta_file <- "matched_sequences_renamed.fasta"
seqs <- readDNAStringSet(fasta_file)
combined_seq <- paste(as.character(seqs), collapse = "")
combined_seq <- DNAString(combined_seq)

# 4. 検出処理
non_cutting <- c()
cat("🔬 スキャン中...\n")

for (enzyme in names(restriction_enzymes)) {
  pattern <- restriction_enzymes[[enzyme]]
  expanded <- expand_iupac(pattern)
  found <- FALSE
  for (p in expanded) {
    if (length(matchPattern(p, combined_seq)) > 0) {
      found <- TRUE
      break
    }
  }
  if (!found) {
    non_cutting <- c(non_cutting, enzyme)
  }
}

# 5. 結果出力
cat("🧬 認識サイトが1箇所もない制限酵素:\n")
print(non_cutting)
