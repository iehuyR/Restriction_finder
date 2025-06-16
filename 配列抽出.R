setwd("~/20250610_newtarget")

library(tidyverse)
library(Biostrings)

# ファイルパスの指定
blast_file <- "20250611_blastp_best_hit_filtered_by_rev_XXX.txt"
fasta_file <- "0124oryctesTotal_nucl.fasta"
output_fasta <- "matched_sequences_renamed.fasta"

# 1. BLAST結果から必要な列を取得
blast_data <- read_tsv(blast_file, show_col_types = FALSE)
sseqid_to_qseqid <- blast_data %>%
  distinct(sseqid, qseqid) %>%
  deframe()  # sseqidを名前、qseqidを値にするリストに変換

# 2. FASTAの読み込み
fasta_seqs <- readDNAStringSet(fasta_file)

# 3. 元ヘッダーの最初の単語だけ取り出す（sseqid部分）
original_ids <- names(fasta_seqs) %>%
  str_split_fixed(" ", 2) %>%
  .[, 1]

# 4. sseqidが一致する配列のみ抽出
matched_indices <- which(original_ids %in% names(sseqid_to_qseqid))
matched_seqs <- fasta_seqs[matched_indices]

# 5. 新しいヘッダーに置き換え（例：seq45474.p1_TC000069-PA）
new_names <- original_ids[matched_indices] %>%
  map_chr(~ paste0(.x, "_", sseqid_to_qseqid[[.x]]))

# 6. 名前を変更して保存
names(matched_seqs) <- new_names
writeXStringSet(matched_seqs, filepath = output_fasta)

cat("✅ ", length(matched_seqs), "件の配列を", output_fasta, "に保存しました。\n")
