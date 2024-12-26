# 创建14个基因文件夹
gene_names <- c("BaiA1", "BaiA2", "BaiB", "BaiCD", "BaiE", "BaiF", "BaiG", "BaiH", "BaiI", "BaiJ", "BaiK", "BaiN", "BaiO", "BaiP")

# 获取当前工作目录
working_dir <- getwd()

# 创建每个基因对应的文件夹
for (gene in gene_names) {
  dir.create(gene, showWarnings = FALSE)  # 如果文件夹已经存在，不会产生警告
}

# 获取工作目录中的所有文件
files <- list.files(working_dir)

# 遍历所有文件并根据基因名进行分类
for (file in files) {
  for (gene in gene_names) {
    # 检查文件名中是否包含基因名
    if (grepl(gene, file)) {
      # 移动文件到对应的文件夹
      file.rename(file, file.path(gene, file))
      cat(paste("Moved file", file, "to folder", gene, "\n"))
      break  # 一旦找到匹配的文件夹，跳出循环
    }
  }
}
