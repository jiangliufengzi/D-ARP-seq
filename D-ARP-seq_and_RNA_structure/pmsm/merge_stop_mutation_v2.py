#!/usr/bin/env python3
"""
================================================================================
脚本名称: merge_stop_mutation_v2.py
功能: 合并 stop 文件和 mutation 文件（pandas 性能优化版）
================================================================================

背景/用途:
  本脚本用于将 call_stop.py 产出的 RT-stop bedgraph 文件与 mutation_FR 产出的
  突变检测文件按位点合并，计算突变率和 RNA 修饰 Signal 值，供下游位点注释和
  修饰信号分析使用。

  典型输入来源:
    - call_stop.py → *.bedgraph (stop 信号)
    - mutation_FR_v4/v5.py → *-mutations.txt (PMSM 突变统计)

  输出用于下游分析:
    - {prefix}-merged.txt 合并表 (含 stop_value + mutation + signal_value)
    - annotate_mRNA_site_ultimate_v3.py 等注释工具输入
    - 修饰位点筛选与 Signal 阈值过滤

  合并键: (chr, position, strand)
    - stop bedgraph 的 end 坐标 (0-based) = mutation 的 position (1-based)
    - 外连接 (outer merge): 保留仅 stop 或仅 mutation 的位点

  派生指标:
    - mutation_rate = (effective_depth - ref_count) / effective_depth
    - signal_value  = (PM + SM + S) / Depth × 100%

  性能优化 (v2):
    - pandas 直接读取 + 向量化计算
    - to_csv 批量写入
    - --large-file-mode 分批处理超大 mutation 文件

================================================================================

工作流程 (按执行顺序):

    ┌─────────────────────────────────────┐
    │  输入: --stop-dir + --mutation-dir  │
    │  可选: --fai-file / --large-file-mode│
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  main() [主入口]                    │
    │  ───────────────────────────────   │
    │  • 解析命令行参数 (argparse)        │
    │  • init_logging() 初始化日志        │
    │  • read_fai_lengths() 染色体长度    │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  find_paired_files()                │
    │  按样本前缀匹配 stop/mutation 对    │
    │  ───────────────────────────────   │
    │  • *.bedgraph ↔ *-mutations.txt     │
    │  • 兼容多种 bedgraph 后缀命名       │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  逐对处理 (循环)                    │
    │  ───────────────────────────────   │
    │  常规模式:                          │
    │  │ read_bedgraph_file_optimized()   │
    │  │ read_mutation_file_optimized()  │
    │  └→ merge_stop_mutation_optimized()│
    │  大文件模式:                        │
    │  └→ process_large_file_in_batches()│
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  finalize_merged_dataframe()        │
    │  ───────────────────────────────   │
    │  • 填充缺失值                       │
    │  • calculate_mutation_rate_vectorized│
    │  • 计算 signal_value                │
    │  • 添加 chr_length                  │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  save_merged_file_optimized()       │
    │  输出 {prefix}-merged.txt           │
    │  日志: logs/merge_stop_mutation_*.log│
    └─────────────────────────────────────┘

================================================================================
核心函数说明 (按执行顺序):

  主流程函数:
  1. main()                           - 主入口，协调整体合并流程
  2. init_logging() / log_message()   - 日志系统
  3. find_paired_files()              - 按样本前缀匹配文件对
  4. read_fai_lengths()               - 读取 .fai/.txt 染色体长度

  数据读取:
  5. read_bedgraph_file_optimized()   - 读取 stop bedgraph (5列)
  6. read_mutation_file_optimized()   - 读取 mutation TSV (14列含 PMSM)

  合并与计算:
  7. merge_stop_mutation_optimized()  - [核心] outer merge + 派生列
     └→ finalize_merged_dataframe()   - 填充缺失值、计算指标、排序
        └→ calculate_mutation_rate_vectorized() - 向量化突变率

  大文件模式:
  8. process_large_file_in_batches()  - 分批读取 mutation，补充 stop-only 位点

  输出:
  9. save_merged_file_optimized()     - to_csv 批量写入合并结果

================================================================================
核心算法说明 - merge_stop_mutation_optimized():

  合并与指标计算:
  ─────────────────────────────────────────────────────────────────────
  步骤    操作                              说明
  ─────────────────────────────────────────────────────────────────────
  1       stop.end → position (1-base)    bedgraph 结束坐标作为合并键
  2       outer merge on (chr,position,   保留 stop-only / mutation-only
          strand)                          和两者都有的位点
  3       fillna 缺失值                   stop/mutation 单侧缺失填 0 或 N
  4       mutation_rate 计算              (depth-N-ref_count)/(depth-N)
  5       signal_value 计算               (PM+SM+S)/depth × 100%
  6       chr_length 映射                 来自 --fai-file (可选)
  7       sort by chr, position, strand   输出排序
  ─────────────────────────────────────────────────────────────────────

  PMSM 值含义 (mutation 文件):
    P  = Pass (正常通过且匹配)       PM = Pass+Mutation
    S  = Stop (RT停止且匹配)        SM = Stop+Mutation

  文件配对规则:
    stop:     {prefix}[_shifted_combined_sorted].bedgraph
    mutation: {prefix}-mutations.txt
    输出:     {prefix}-merged.txt

================================================================================
输入格式:

  Stop 文件 (*.bedgraph):
  ─────────────────────────────────────────────────────────────────────
  列          坐标系    说明
  ─────────────────────────────────────────────────────────────────────
  chr         -         序列名
  start       0-based   起始坐标
  end         0-based   结束坐标 (= 合并用的 1-base position)
  stop_value  -         RT-stop 信号值
  strand      -         链方向 (+/-)
  ─────────────────────────────────────────────────────────────────────

  Mutation 文件 (*-mutations.txt):
  ─────────────────────────────────────────────────────────────────────
  列          说明
  ─────────────────────────────────────────────────────────────────────
  chr         序列名
  position    1-based 位置
  ref_base    参考碱基
  strand      链方向
  depth       覆盖深度
  A/C/G/T/N   碱基计数
  P/PM/S/SM   PMSM 统计值
  ─────────────────────────────────────────────────────────────────────

  染色体长度 (--fai-file, 可选):
    - .fai 标准格式 (取前 2 列) 或 .txt 两列 (seq_name, length)

================================================================================
输出格式:

  合并结果 ({prefix}-merged.txt):
  ─────────────────────────────────────────────────────────────────────
  列名            类型    说明
  ─────────────────────────────────────────────────────────────────────
  chr             str     序列名
  position        int     1-based 位置
  ref_base        str     参考碱基
  strand          str     链方向
  stop_value      float   RT-stop 信号 (无 stop 时为 0)
  depth           int     覆盖深度
  A/C/G/T/N       int     碱基计数
  P/PM/S/SM       int     PMSM 统计
  mutation_rate   float   突变率 (4位小数)
  signal_value    float   修饰信号 % (2位小数)
  chr_length      int     染色体长度 (无 fai 时为 0)
  ─────────────────────────────────────────────────────────────────────

  日志: {output_dir}/logs/merge_stop_mutation_{timestamp}.log

================================================================================
命令行参数:

  必需/主要参数:
    --stop-dir              stop bedgraph 目录 (默认: 8.stop)
    --mutation-dir          mutation 文件目录 (默认: 9.mutation)
    --output-dir            输出目录 (默认: 10.merged)

  可选参数:
    --fai-file              染色体长度文件 (.fai 或 .txt)
    --large-file-mode       启用大文件分批处理模式
    --batch-size            分批大小 (默认: 1000000)

================================================================================
使用示例:

  1. 基本合并:
      python /home/pf/14T/scripts/merge_stop_mutation_v2.py \\
         --stop-dir 8.stop --mutation-dir 9.mutation --output-dir 10.merged

  2. 附带染色体长度:
      python /home/pf/14T/scripts/merge_stop_mutation_v2.py \\
         --stop-dir 8.stop --mutation-dir 9.mutation \\
         --output-dir 10.merged \\
         --fai-file reference.fa.fai

  3. 大文件分批模式:
      python /home/pf/14T/scripts/merge_stop_mutation_v2.py \\
         --stop-dir 8.stop --mutation-dir 9.mutation \\
         --output-dir 10.merged \\
         --large-file-mode --batch-size 500000

  4. 完整参数示例:
      python /home/pf/14T/scripts/merge_stop_mutation_v2.py \\
         --stop-dir /data/8.stop \\
         --mutation-dir /data/9.mutation \\
         --output-dir /data/10.merged \\
         --fai-file /ref/genome.fa.fai \\
         --large-file-mode --batch-size 1000000

================================================================================
依赖:
  - Python 3.6+
  - pandas, numpy
  - 标准库: os, glob, argparse, logging, pathlib

================================================================================
注意事项:
  1. v3 (Polars 版) 输出列兼容 v2，大数据集建议用 merge_stop_mutation_v3.py
  2. 合并使用 outer join，单侧缺失位点会保留并用默认值填充
  3. bedgraph 的 end 坐标必须与 mutation position 在同一坐标系下对应
  4. mutation 文件缺少 PMSM 列时自动补 0
  5. --large-file-mode 先分批处理 mutation，最后补充 stop-only 位点

================================================================================
"""

import os
import pandas as pd
import numpy as np
import glob
from collections import defaultdict
import argparse
import warnings
import logging
import datetime
from pathlib import Path
warnings.filterwarnings('ignore')

# 全局日志文件路径
_log_path = None


def init_logging(output_dir):
    """初始化日志系统"""
    global _log_path
    
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(exist_ok=True, parents=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    _log_path = log_dir / f"merge_stop_mutation_{timestamp}.log"


def log_message(message):
    """记录日志消息到控制台和文件"""
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_entry = f"[{timestamp}] {message}"
    
    # 打印到控制台
    print(log_entry)
    
    # 写入日志文件
    if _log_path:
        try:
            with open(_log_path, 'a', encoding='utf-8') as f:
                f.write(log_entry + '\n')
        except Exception as e:
            print(f"警告: 无法写入日志文件: {e}")


def read_fai_lengths(fai_file):
    """
    读取染色体长度文件，根据文件扩展名自动判断格式
    - .fai文件: 标准fai格式，使用前2列（序列名、长度）
    - .txt文件: 简化格式，2列（序列名、长度）  
    - 其他扩展名: 自动检测列数
    """
    import os
    
    try:
        file_ext = os.path.splitext(fai_file)[1].lower()
        
        if file_ext == '.fai':
            # 标准fai格式: 至少5列，使用前2列
            log_message(f"  检测到.fai格式文件，使用标准fai格式读取")
            fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                names=['seq_name', 'length', 'offset', 'linebases', 'linewidth'],
                                usecols=[0, 1])
        elif file_ext == '.txt':
            # 简化txt格式: 2列
            log_message(f"  检测到.txt格式文件，使用2列格式读取")
            fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                names=['seq_name', 'length'])
        else:
            # 其他扩展名：自动检测列数
            log_message(f"  未知扩展名{file_ext}，自动检测列数...")
            with open(fai_file, 'r') as f:
                first_line = f.readline().strip()
                n_cols = len(first_line.split('\t'))
            
            if n_cols == 2:
                log_message(f"  检测到{n_cols}列，使用2列格式读取")
                fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                    names=['seq_name', 'length'])
            elif n_cols >= 5:
                log_message("使用标准fai格式读取")
                fai_df = pd.read_csv(fai_file, sep='\t', header=None, 
                                    names=['seq_name', 'length', 'offset', 'linebases', 'linewidth'],
                                    usecols=[0, 1])
            else:
                raise ValueError(f"不支持的文件格式，列数: {n_cols}")
            
        return dict(zip(fai_df['seq_name'], fai_df['length']))
    except Exception as e:
        log_message(f"警告: 读取染色体长度文件时出错: {e}")
        return {}


def find_paired_files(stop_dir, mutation_dir):
    """
    找到配对的stop和mutation文件
    """
    stop_files = glob.glob(os.path.join(stop_dir, "*.bedgraph"))
    mutation_files = glob.glob(os.path.join(mutation_dir, "*-mutations.txt"))

    def normalize_stop_prefix(file_path):
        basename = os.path.basename(file_path)
        known_suffixes = [
            "_shifted_combined_sorted.bedgraph",
            "_shifted_combined.bedgraph",
            ".bedgraph",
        ]
        for suffix in known_suffixes:
            if basename.endswith(suffix):
                return basename[:-len(suffix)]
        return os.path.splitext(basename)[0]

    def normalize_mutation_prefix(file_path):
        basename = os.path.basename(file_path)
        suffix = "-mutations.txt"
        if basename.endswith(suffix):
            return basename[:-len(suffix)]
        return os.path.splitext(basename)[0]
    
    # 提取样本前缀，兼容call_stop.py和mutation_FR_v4.py的默认输出文件名
    stop_prefixes = {}
    for file in stop_files:
        prefix = normalize_stop_prefix(file)
        stop_prefixes[prefix] = file
    
    mutation_prefixes = {}
    for file in mutation_files:
        prefix = normalize_mutation_prefix(file)
        mutation_prefixes[prefix] = file
    
    # 找到配对的文件
    paired_files = []
    for prefix in stop_prefixes:
        if prefix in mutation_prefixes:
            paired_files.append((prefix, stop_prefixes[prefix], mutation_prefixes[prefix]))
    
    return paired_files


def read_bedgraph_file_optimized(file_path):
    """
    优化版：使用pandas直接读取bedgraph文件
    """
    try:
        df = pd.read_csv(file_path, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'stop_value', 'strand'],
                        dtype={'chr': str, 'start': int, 'end': int, 
                              'stop_value': float, 'strand': str})
        df['position_1base'] = df['end']
        return df
    except Exception as e:
        log_message(f"  警告: 读取bedgraph文件时出错，使用备用方法: {e}")
        # 备用方法：处理可能的空文件或格式问题
        return pd.DataFrame(columns=['chr', 'start', 'end', 'stop_value', 'strand', 'position_1base'])


def read_mutation_file_optimized(file_path):
    """
    优化版：使用pandas直接读取mutation文件
    修复：使用更宽松的数据类型，避免读取失败
    """
    try:
        # 先检查文件有多少列
        with open(file_path, 'r') as f:
            header_line = f.readline()
            sample_line = f.readline()
            if sample_line:
                n_cols = len(sample_line.strip().split('\t'))
            else:
                n_cols = 14  # 默认列数
        
        log_message(f"  检测到mutation文件有 {n_cols} 列")
        
        # 根据列数选择合适的列名，但不强制dtype（避免类型错误）
        if n_cols >= 14:
            col_names = ['chr', 'position_1base', 'ref_base', 'strand', 'depth', 
                        'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM']
        else:
            col_names = ['chr', 'position_1base', 'ref_base', 'strand', 'depth', 
                        'A', 'C', 'G', 'T', 'N']
        
        # 读取文件 - 不指定dtype，让pandas自动推断
        log_message(f"  正在读取mutation文件...")
        read_cols = min(n_cols, len(col_names))
        df = pd.read_csv(file_path, sep='\t', skiprows=1,
                        names=col_names[:read_cols],
                        usecols=list(range(read_cols)))
        
        log_message(f"  成功读取 {len(df):,} 行mutation数据")
        
        # 如果缺少P, PM, S, SM列，添加默认值
        for col in ['P', 'PM', 'S', 'SM']:
            if col not in df.columns:
                df[col] = 0
        
        # 转换数据类型（容错处理）
        try:
            # 数值列转换为数值类型
            numeric_cols = ['position_1base', 'depth', 'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM']
            for col in numeric_cols:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
        except Exception as convert_error:
            log_message(f"  警告: 数据类型转换时出现问题: {convert_error}")
        
        return df
    except Exception as e:
        log_message(f"  ❌ 错误: 读取mutation文件失败: {e}")
        log_message(f"  返回空DataFrame，这会导致所有mutation数据为默认值！")
        return pd.DataFrame(columns=['chr', 'position_1base', 'ref_base', 'strand', 'depth',
                                    'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM'])


def calculate_mutation_rate_vectorized(df):
    """
    向量化计算突变率（性能提升5-10倍）
    """
    # 创建ref_count列 - 使用numpy.select进行向量化条件选择
    conditions = [
        df['ref_base'] == 'A',
        df['ref_base'] == 'C',
        df['ref_base'] == 'G',
        df['ref_base'] == 'T'
    ]
    choices = [df['A'], df['C'], df['G'], df['T']]
    ref_count = np.select(conditions, choices, default=0)
    
    # 计算有效深度
    effective_depth = df['depth'] - df['N']
    
    # 创建掩码：有效数据的条件
    total_bases = df['A'] + df['C'] + df['G'] + df['T']
    valid_mask = (df['depth'] > 0) & (effective_depth > 0) & (total_bases > 0)
    
    # 初始化突变率列
    df['mutation_rate'] = 0.0
    
    # 向量化计算突变率
    df.loc[valid_mask, 'mutation_rate'] = (
        (effective_depth[valid_mask] - ref_count[valid_mask]) / effective_depth[valid_mask]
    )
    
    # 确保突变率在0-1之间
    df['mutation_rate'] = df['mutation_rate'].clip(lower=0, upper=1)
    
    return df


def finalize_merged_dataframe(merged_df, chromosome_lengths=None):
    """
    填充缺失值并计算输出所需的派生列。
    """
    # 批量填充缺失值（使用更高效的方式）
    log_message("  正在填充缺失值...")
    # 使用字典批量填充，性能更好
    fill_values = {
        'stop_value': 0,
        'ref_base': 'N',
        'depth': 0,
        'A': 0,
        'C': 0,
        'G': 0,
        'T': 0,
        'N': 0,
        'P': 0,
        'PM': 0,
        'S': 0,
        'SM': 0
    }
    for col, value in fill_values.items():
        if col not in merged_df.columns:
            merged_df[col] = value
    merged_df = merged_df.fillna(fill_values)

    if 'stop_value' in merged_df.columns:
        merged_df['stop_value'] = pd.to_numeric(merged_df['stop_value'], errors='coerce').fillna(0)
    
    # 确保整数列的类型正确
    int_cols = ['depth', 'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM', 'position']
    for col in int_cols:
        if col in merged_df.columns:
            merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').fillna(0).astype(np.int32)
    
    log_message("  正在计算突变率（向量化）...")
    merged_df = calculate_mutation_rate_vectorized(merged_df)
    
    log_message("  正在计算Signal值（向量化）...")
    # 向量化计算Signal值: (PM + SM + S) / Depth × 100%
    signal_sum = merged_df['PM'] + merged_df['SM'] + merged_df['S']
    depth_mask = merged_df['depth'] > 0
    merged_df['signal_value'] = 0.0
    merged_df.loc[depth_mask, 'signal_value'] = (
        signal_sum[depth_mask] / merged_df.loc[depth_mask, 'depth'] * 100.0
    )
    
    # 添加染色体长度（向量化）
    if chromosome_lengths:
        log_message("  正在添加染色体长度信息...")
        merged_df['chr_length'] = merged_df['chr'].map(chromosome_lengths).fillna(0).astype(np.int32)
    else:
        merged_df['chr_length'] = 0
    
    # 使用更高效的排序
    log_message("  正在排序...")
    merged_df = merged_df.sort_values(['chr', 'position', 'strand'], ignore_index=True)
    
    return merged_df


def merge_stop_mutation_optimized(stop_df, mutation_df, chromosome_lengths=None):
    """
    优化版合并函数
    """
    # 重命名列以便合并
    stop_df_renamed = stop_df.rename(columns={'position_1base': 'position'})
    mutation_df_renamed = mutation_df.rename(columns={'position_1base': 'position'})
    
    # 使用pandas高效外连接合并
    log_message("  正在合并数据...")
    merged_df = pd.merge(
        mutation_df_renamed, 
        stop_df_renamed[['chr', 'position', 'strand', 'stop_value']], 
        on=['chr', 'position', 'strand'], 
        how='outer'
    )
    
    return finalize_merged_dataframe(merged_df, chromosome_lengths)


def save_merged_file_optimized(merged_df, output_path):
    """
    优化版：使用to_csv批量写入（速度提升10倍以上）
    """
    # 确保列顺序正确
    columns = ['chr', 'position', 'ref_base', 'strand', 'stop_value', 
               'depth', 'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM',
               'mutation_rate', 'signal_value', 'chr_length']
    
    # 格式化数值列
    output_df = merged_df[columns].copy()
    output_df['mutation_rate'] = output_df['mutation_rate'].round(4)
    output_df['signal_value'] = output_df['signal_value'].round(2)
    output_df['stop_value'] = output_df['stop_value'].astype(float)
    
    # 一次性写入文件
    output_df.to_csv(output_path, sep='\t', index=False, float_format='%.4g')


def process_large_file_in_batches(stop_file, mutation_file, chromosome_lengths, batch_size=1000000):
    """
    对于超大文件，分批处理以减少内存使用
    """
    # 读取stop文件（通常较小）
    stop_df = read_bedgraph_file_optimized(stop_file)
    stop_df_renamed = stop_df.rename(columns={'position_1base': 'position'})
    stop_lookup = stop_df_renamed[['chr', 'position', 'strand', 'stop_value']]
    
    # 分批处理mutation文件
    merged_chunks = []
    mutation_key_chunks = []
    
    try:
        # 获取文件总行数（用于进度显示）
        with open(mutation_file, 'r') as f:
            header_line = f.readline()
            sample_line = f.readline()
            n_cols = len(sample_line.strip().split('\t')) if sample_line else 14
            total_lines = 1 + sum(1 for line in f) if sample_line else 0
        
        n_batches = (total_lines + batch_size - 1) // batch_size
        col_names = ['chr', 'position_1base', 'ref_base', 'strand', 'depth',
                     'A', 'C', 'G', 'T', 'N', 'P', 'PM', 'S', 'SM']
        read_cols = min(n_cols, len(col_names))
        
        # 分批读取和处理
        for i, chunk in enumerate(pd.read_csv(mutation_file, sep='\t', skiprows=1, 
                                             chunksize=batch_size, header=None,
                                             names=col_names[:read_cols],
                                             usecols=list(range(read_cols)))):
            log_message(f"    处理批次 {i+1}/{n_batches}...")
            
            for col in ['P', 'PM', 'S', 'SM']:
                if col not in chunk.columns:
                    chunk[col] = 0
            
            chunk_renamed = chunk.rename(columns={'position_1base': 'position'})
            mutation_key_chunks.append(chunk_renamed[['chr', 'position', 'strand']].copy())

            # 当前批次只做left merge，最后再补充没有mutation记录的stop-only位点。
            merged_chunk = pd.merge(
                chunk_renamed,
                stop_lookup,
                on=['chr', 'position', 'strand'],
                how='left'
            )
            merged_chunk = finalize_merged_dataframe(merged_chunk, chromosome_lengths)
            merged_chunks.append(merged_chunk)
    
    except Exception as e:
        log_message(f"    批处理出错，尝试常规方法: {e}")
        mutation_df = read_mutation_file_optimized(mutation_file)
        return merge_stop_mutation_optimized(stop_df, mutation_df, chromosome_lengths)

    if mutation_key_chunks:
        all_mutation_keys = pd.concat(mutation_key_chunks, ignore_index=True).drop_duplicates()
        stop_only = pd.merge(
            stop_lookup,
            all_mutation_keys,
            on=['chr', 'position', 'strand'],
            how='left',
            indicator=True
        )
        stop_only = stop_only[stop_only['_merge'] == 'left_only'].drop(columns=['_merge'])
    else:
        stop_only = stop_lookup

    if not stop_only.empty:
        log_message("    添加仅存在于stop文件中的位点...")
        merged_chunks.append(finalize_merged_dataframe(stop_only, chromosome_lengths))
    
    # 合并所有批次
    if merged_chunks:
        final_df = pd.concat(merged_chunks, ignore_index=True)
        final_df = final_df.sort_values(['chr', 'position', 'strand'], ignore_index=True)
        return final_df
    else:
        return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(description='合并stop和mutation文件（性能优化版）')
    parser.add_argument('--stop-dir', default='8.stop', help='stop文件夹路径')
    parser.add_argument('--mutation-dir', default='9.mutation', help='mutation文件夹路径')
    parser.add_argument('--output-dir', default='10.merged', help='输出文件夹路径')
    parser.add_argument('--fai-file', required=False, help='染色体fasta索引文件(.fai)路径，用于获取染色体长度信息')
    parser.add_argument('--large-file-mode', action='store_true', help='启用大文件模式（分批处理）')
    parser.add_argument('--batch-size', type=int, default=1000000, help='大文件模式下的批次大小')
    
    args = parser.parse_args()
    
    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 初始化日志系统
    init_logging(args.output_dir)
    
    # 读取染色体长度信息
    chromosome_lengths = None
    if args.fai_file:
        try:
            chromosome_lengths = read_fai_lengths(args.fai_file)
            log_message(f"从fai文件读取到 {len(chromosome_lengths)} 个染色体长度信息")
        except Exception as e:
            log_message(f"警告: 无法读取fai文件 {args.fai_file}: {str(e)}")
            log_message("将继续处理，但不包含染色体长度信息")
    
    # 找到配对的文件
    paired_files = find_paired_files(args.stop_dir, args.mutation_dir)
    
    log_message(f"找到 {len(paired_files)} 对配对文件")
    if args.large_file_mode:
        log_message(f"大文件模式已启用，批次大小: {args.batch_size}")
    log_message("")
    
    # 记录总体统计
    total_processed = 0
    total_errors = 0
    
    for idx, (prefix, stop_file, mutation_file) in enumerate(paired_files, 1):
        log_message(f"[{idx}/{len(paired_files)}] 处理: {prefix}")
        log_message(f"  Stop文件: {os.path.basename(stop_file)}")
        log_message(f"  Mutation文件: {os.path.basename(mutation_file)}")
        
        try:
            output_file = os.path.join(args.output_dir, f"{prefix}-merged.txt")
            if args.large_file_mode:
                # 复用v3的SQLite磁盘后端；旧版“分块”最终仍会concat全部块。
                from merge_stop_mutation_v3 import process_large_file
                process_large_file(
                    stop_file, mutation_file, chromosome_lengths, output_file,
                    args.batch_size
                )
                log_message(f"  ✓ 输出: {output_file}")
                total_processed += 1
                log_message("")
                continue
            else:
                # 常规模式：一次性加载
                stop_df = read_bedgraph_file_optimized(stop_file)
                mutation_df = read_mutation_file_optimized(mutation_file)
                
                log_message(f"  Stop数据: {len(stop_df):,} 行")
                log_message(f"  Mutation数据: {len(mutation_df):,} 行")
                
                # 合并数据
                merged_df = merge_stop_mutation_optimized(stop_df, mutation_df, chromosome_lengths)
            
            log_message(f"  合并后: {len(merged_df):,} 行")
            
            # 保存结果
            save_merged_file_optimized(merged_df, output_file)
            log_message(f"  ✓ 输出: {output_file}")
            
            total_processed += 1
            
        except Exception as e:
            log_message(f"  ✗ 错误: {str(e)}")
            total_errors += 1
        
        log_message("")  # 空行分隔
    
    log_message("=" * 50)
    log_message(f"处理完成!")
    log_message(f"成功: {total_processed} 个文件")
    if total_errors > 0:
        log_message(f"失败: {total_errors} 个文件")
    log_message(f"结果保存在 {args.output_dir} 文件夹中")
    log_message(f"日志文件: {_log_path}")


if __name__ == "__main__":
    main()
