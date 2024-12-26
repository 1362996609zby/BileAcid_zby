# 加载所需的库
if (!require("openxlsx")) install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# 获取当前目录下所有 .txt 文件
txt_files <- list.files(path = ".", pattern = "\\.txt$", full.names = TRUE)

# 获取当前目录名
current_dir_name <- basename(getwd())

# 初始化存储数据的列表
result_data <- data.frame(
  `.txt列名` = character(),
  stringsAsFactors = FALSE
)

# 遍历所有 .txt 文件
data_counts <- c()  # 用于存储每个文件的数据个数
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
    
    # 统计数据行数
    data_counts <- c(data_counts, length(data_lines))
    result_data <- rbind(result_data, data.frame(`.txt列名` = basename(file), stringsAsFactors = FALSE))
  } else {
    cat("文件无目标表头:", basename(file), "\n")
    data_counts <- c(data_counts, 0)  # 无有效数据计数为0
    result_data <- rbind(result_data, data.frame(`.txt列名` = basename(file), stringsAsFactors = FALSE))
  }
}

# 添加当前目录名列
result_data[[current_dir_name]] <- data_counts

# 将结果写入一个 Excel 文件
output_file <- "file_data_summary_reformatted.xlsx"
write.xlsx(result_data, file = output_file, rowNames = FALSE)
cat("统计结果已保存到:", output_file, "\n")
