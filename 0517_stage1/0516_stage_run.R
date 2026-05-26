# ============================================================
# 主控脚本：依次运行多个数值模拟代码
# 说明：
# 1. 每个子脚本在独立环境中运行，避免 rm(list = ls()) 影响主控脚本
# 2. 每个子脚本的 cat()、print()、message() 等输出会正常显示
# 3. 某个脚本出错后，不会中断后续脚本运行
# ============================================================

rm(list = ls())

# ----------------------------
# 1. 基本设置
# ----------------------------

start_time <- Sys.time()
main_wd <- getwd()

cat("\n============================================================\n")
cat("数值模拟开始运行\n")
cat("开始时间：", format(start_time), "\n")
cat("当前工作目录：", main_wd, "\n")
cat("============================================================\n\n")

# 如有需要，可手动设置代码所在文件夹
# setwd("D:/你的代码所在文件夹")
# main_wd <- getwd()

script_files <- c(
# "0516_stage1_single_1.R",
# "0516_stage1_single_2.R",
#"0517_stage1_double.R"
   "0518_stage2_1.R",
   "0518_stage2_2.R"
)

# 用于记录每个脚本的运行结果
run_summary <- data.frame(
  文件名 = character(),
  状态 = character(),
  用时_分钟 = numeric(),
  错误信息 = character(),
  stringsAsFactors = FALSE
)

# ----------------------------
# 2. 依次运行代码文件
# ----------------------------

for (i in seq_along(script_files)) {
  
  current_script <- script_files[i]
  script_path <- file.path(main_wd, current_script)
  
  cat("\n------------------------------------------------------------\n")
  cat("开始运行第", i, "个文件：", current_script, "\n")
  cat("------------------------------------------------------------\n\n")
  
  script_start <- Sys.time()
  status <- "成功"
  error_msg <- ""
  
  tryCatch(
    {
      if (!file.exists(script_path)) {
        stop("找不到文件：", script_path)
      }
      
      # 每个子脚本使用独立环境运行
      # 这样即使子脚本中有 rm(list = ls())，也不会清空主控脚本变量
      script_env <- new.env(parent = globalenv())
      
      source(
        file = script_path,
        local = script_env,
        encoding = "UTF-8",
        echo = FALSE,
        print.eval = TRUE
      )
    },
    error = function(e) {
      status <<- "失败"
      error_msg <<- conditionMessage(e)
      
      cat("\n运行出错：", current_script, "\n")
      cat("错误信息：", error_msg, "\n")
    },
    finally = {
      # 如果子脚本中修改了工作目录，这里恢复到主控脚本所在目录
      setwd(main_wd)
    }
  )
  
  script_end <- Sys.time()
  used_time <- round(difftime(script_end, script_start, units = "mins"), 2)
  
  if (status == "成功") {
    cat("\n运行完成：", current_script, "\n")
    cat("用时：", used_time, "分钟\n")
  } else {
    cat("\n运行失败：", current_script, "\n")
    cat("已继续运行后续脚本。\n")
    cat("用时：", used_time, "分钟\n")
  }
  
  run_summary <- rbind(
    run_summary,
    data.frame(
      文件名 = current_script,
      状态 = status,
      用时_分钟 = as.numeric(used_time),
      错误信息 = error_msg,
      stringsAsFactors = FALSE
    )
  )
}

# ----------------------------
# 3. 结束信息
# ----------------------------

end_time <- Sys.time()

cat("\n============================================================\n")
cat("全部代码运行结束\n")
cat("结束时间：", format(end_time), "\n")
cat("总用时：",
    round(difftime(end_time, start_time, units = "mins"), 2),
    "分钟\n")
cat("============================================================\n\n")

cat("各文件运行汇总：\n")
print(run_summary)