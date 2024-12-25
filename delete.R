# 获取当前目录下所有 .txt 文件
txt_files <- list.files(path = ".", pattern = "\\.txt$", full.names = TRUE)

# 遍历每个文件
for (file in txt_files) {
  # 提取文件名和路径
  file_path <- dirname(file)
  file_name <- basename(file)
  
  # 删除 '_translated_hmmsearch' 后缀
  new_name <- sub("_translated_hmmsearch", "", file_name)
  
  # 构建新文件路径
  new_file_path <- file.path(file_path, new_name)
  
  # 重命名文件
  file.rename(file, new_file_path)
  cat("文件已重命名:", file_name, "->", new_name, "\n")
}
