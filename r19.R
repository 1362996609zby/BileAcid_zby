# 加载所需的库
if (!require("openxlsx")) install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# 获取当前目录下所有 .txt 文件
txt_files <- list.files(path = ".", pattern = "\\.txt$", full.names = TRUE)

# 初始化存储数据的列表
combined_data <- NULL

# 定义标准列名
default_colnames <- c(
  "target_name", "accession", "query_name", "query_accession", "E-value_full", 
  "score_full", "bias_full", "E-value_domain", "score_domain", "bias_domain", 
  "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description"
)

# 遍历所有 .txt 文件
for (file in txt_files) {
  # 读取文件内容
  cat("正在处理文件:", basename(file), "\n")
  file_content <- readLines(file)
  
  # 检查文件是否包含目标表头
  header_index <- grep("^# target name.*description of target$", file_content)
  if (length(header_index) > 0) {
    # 提取表头及其下方的数据
    data_start <- header_index + 2 # 数据从表头下两行开始
    data_lines <- file_content[data_start:length(file_content)]
    
    # 移除空行和注释行
    data_lines <- data_lines[!grepl("^#|^\\s*$", data_lines)]
    
    # 如果存在数据，转换为数据框
    if (length(data_lines) > 0) {
      # 分割行内容为列，使用空格作为分隔符（忽略空格数量）
      data <- do.call(rbind, strsplit(trimws(data_lines), "\\s+"))
      
      # 转换为数据框
      data <- as.data.frame(data, stringsAsFactors = FALSE)
      
      # 修正列数到19
      if (ncol(data) > 19) {
        data[, 19] <- apply(data[, 19:ncol(data)], 1, paste, collapse = " ")
        data <- data[, 1:19]
      } else if (ncol(data) < 19) {
        data <- cbind(data, matrix(NA, nrow = nrow(data), ncol = 19 - ncol(data)))
      }
      
      # 设置列名
      colnames(data) <- default_colnames
      
      # 添加文件来源列
      data$file <- basename(file)  # 添加来源文件列
      
      # 合并到总数据中
      combined_data <- rbind(combined_data, data)
      cat("文件处理完成:", basename(file), "\n")
    } else {
      cat("文件无数据内容:", basename(file), "\n")
    }
  } else {
    cat("文件无目标表头:", basename(file), "\n")
  }
}

# 如果有结果，将其写入一个 Excel 文件
if (!is.null(combined_data)) {
  # 写入 Excel 文件
  write.xlsx(combined_data, file = "hmmsearch_results_combined.xlsx", rowNames = FALSE)
  cat("所有数据已合并并保存到 hmmsearch_results_combined.xlsx\n")
} else {
  cat("没有发现符合条件的文件或数据。\n")
}
