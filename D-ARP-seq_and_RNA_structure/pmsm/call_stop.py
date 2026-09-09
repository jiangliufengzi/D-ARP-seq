#!/usr/bin/env python3
"""
================================================================================
脚本名称: call_stop.py
功能: 从BAM文件中提取reads末端位点信息，用于停顿位点(stop site)分析
================================================================================

背景/用途:
  本脚本用于RNA结构探测(如DMS-seq, SHAPE-seq)、核糖体profiling、翻译组学等
  分析流程中，从比对后的BAM文件提取reads的5'端或3'端位置信息。
  
  核心原理：
  - reads的5'端或3'端对应逆转录酶停顿位点或RNA修饰位点
  - 通过shift操作校正测序接头/引物偏移
  - 分别统计正负链覆盖度，保留链特异性信息

  典型输入来源：
  - 比对软件(Bowtie2/STAR等)产生的排序BAM文件
  - 去重后的BAM文件（如 *.rmdup.bam）

  输出用途：
  - 停顿位点(stop site)分析和可视化
  - mutation-truncation分析的truncation部分
  - RNA修饰位点的鉴定（如m1A, m5C等会导致RT停顿）
  - 核糖体A/P/E位点分析（Ribo-seq）

  依赖工具：
  - bedtools (bamtobed, shift, genomecov)
  - genome文件: 染色体名称和长度的Tab分隔文件 (chrNameLength.txt)

================================================================================
工作流程 (按执行顺序):

    ┌─────────────────────────────────────┐
    │  输入: BAM文件目录                   │
    │  - 包含 *.bam 文件                  │
    │  - 需要genome size文件              │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  步骤1: 参数解析与初始化             │
    │  ───────────────────────────────   │
    │  • 解析命令行参数                   │
    │  • 设置日志 (output_dir/logs/目录) │
    │  • 扫描BAM文件                      │
    │    - 支持通配符: --pattern          │
    │    - 支持前缀匹配: --pattern-begin  │
    │    - 支持中间匹配: --pattern-mid    │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  步骤2: process_bam_parallel()      │
    │  并行/串行处理调度                   │
    │  ───────────────────────────────   │
    │  • workers=1时串行处理              │
    │  • workers>1时使用ThreadPoolExecutor│
    │  │ ┌────────────────────────────┐  │
    │  └→│ process_single_bam() × N   │  │
    │    │ 单个BAM文件完整处理流程    │  │
    │    └────────────────────────────┘  │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────────────┐
    │  process_single_bam() 单文件处理流程         │
    │  ─────────────────────────────────────────  │
    │                                             │
    │  ┌─────────────────────────────────────┐   │
    │  │ Step 1: BAM → BED                   │   │
    │  │ • bedtools bamtobed                 │   │
    │  │ • 输出: {sample}.bed                │   │
    │  └─────────────┬───────────────────────┘   │
    │                │                            │
    │                ↓                            │
    │  ┌─────────────────────────────────────┐   │
    │  │ Step 2: Shift操作                   │   │
    │  │ • bedtools shift                    │   │
    │  │ • 正链偏移: --plus-shift (默认0)    │   │
    │  │ • 负链偏移: --minus-shift (默认0)   │   │
    │  │ • 输出: {sample}_shifted.bed        │   │
    │  └─────────────┬───────────────────────┘   │
    │                │                            │
    │                ↓                            │
    │  ┌─────────────────────────────────────┐   │
    │  │ Step 3: Plus链覆盖度                │   │
    │  │ • bedtools genomecov -strand +      │   │
    │  │ • -5 或 -3 (统计5'/3'端)            │   │
    │  │ • 输出: {sample}_shifted_plus.bg    │   │
    │  └─────────────┬───────────────────────┘   │
    │                │                            │
    │                ↓                            │
    │  ┌─────────────────────────────────────┐   │
    │  │ Step 4: Minus链覆盖度               │   │
    │  │ • bedtools genomecov -strand -      │   │
    │  │ • -5 或 -3 (统计5'/3'端)            │   │
    │  │ • 输出: {sample}_shifted_minus.bg   │   │
    │  └─────────────┬───────────────────────┘   │
    │                │                            │
    │                ↓                            │
    │  ┌─────────────────────────────────────┐   │
    │  │ Step 5: 合并正负链并排序            │   │
    │  │ • 添加strand列 (+/-)                │   │
    │  │ • LC_COLLATE=C sort排序             │   │
    │  │ • 输出: {sample}_shifted_combined_  │   │
    │  │         sorted.bedgraph             │   │
    │  └─────────────────────────────────────┘   │
    │                                             │
    │  (可选) 删除中间文件 (默认删除)             │
    └──────────────┬──────────────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  输出文件汇总                        │
    │  ───────────────────────────────   │
    │  输出目录/                          │
    │  ├── sample1_shifted_combined_     │
    │  │   sorted.bedgraph               │
    │  ├── sample2_shifted_combined_     │
    │  │   sorted.bedgraph               │
    │  └── ...                           │
    │                                     │
    │  输出目录/logs/                     │
    │  └── bedtools_processing_          │
    │      YYYYMMDD_HHMMSS.log           │
    └─────────────────────────────────────┘

================================================================================
核心函数说明 (按执行顺序):

  1. setup_logger()           - 配置日志系统，同时输出到文件和控制台
  2. get_duration()           - 计算耗时，格式化为 HH:MM:SS.mmm
  3. run_command()            - 执行shell命令的封装函数
  4. process_bam_parallel()   - 主入口，扫描文件并调度并行/串行处理
     └→ process_single_bam()  - 单文件处理核心函数 (内部调用)
        ├→ bedtools bamtobed  - BAM转BED
        ├→ bedtools shift     - 坐标偏移校正
        ├→ bedtools genomecov - 计算覆盖度 (正负链分别)
        └→ sort               - 合并排序

================================================================================
输入格式:

  BAM文件要求:
  - 格式: 标准BAM格式 (*.bam)
  - 建议已排序

  Genome文件格式 (chrNameLength.txt):
  - Tab分隔，两列: 染色体名称 \t 长度
  - 示例:
      chr1    248956422
      chr2    242193529
      ...

================================================================================
输出格式:

  最终输出文件 (*_shifted_combined_sorted.bedgraph):
  - 5列Tab分隔: chrom, start, end, count, strand
  - 示例:
      chr1    10000    10001    5    +
      chr1    10050    10051    3    -

  日志文件 (output_dir/logs/bedtools_processing_*.log):
  - 记录每个步骤的开始/完成时间
  - 记录处理进度和错误信息

================================================================================
使用示例:

  # 1. 最小运行 - 处理当前目录下所有BAM文件
    python /home/pf/14T/scripts/call_stop.py -g /path/to/chrNameLength.txt

  # 2. 指定输入输出目录
    python /home/pf/14T/scripts/call_stop.py -i 3.mapped -o 4.stop_sites -g genome.txt

  # 3. 使用文件名前缀筛选特定样本
    python /home/pf/14T/scripts/call_stop.py -i bams -o output -g genome.txt --pattern-begin "treated"

  # 4. 使用文件名中间部分筛选
    python /home/pf/14T/scripts/call_stop.py -i bams -o output -g genome.txt --pattern-mid "rep1"

  # 5. 组合前缀和中间部分筛选
    python /home/pf/14T/scripts/call_stop.py -o results --pattern-begin "sample" --pattern-mid "treated"

  # 6. 并行处理 (4个线程)
    python /home/pf/14T/scripts/call_stop.py -i bams -o output -g genome.txt -w 4

  # 7. 统计3'端而非5'端
    python /home/pf/14T/scripts/call_stop.py -i bams -o output -g genome.txt --end-type 3

  # 8. 自定义shift值 (如Ribo-seq的P-site校正)
    python /home/pf/14T/scripts/call_stop.py -i bams -o output -g genome.txt --plus-shift -12 --minus-shift 12

  # 9. 保留所有中间文件用于调试
    python /home/pf/14T/scripts/call_stop.py -i bams -o output -g genome.txt --keep-temp

  # 10. 完整参数示例
    python /home/pf/14T/scripts/call_stop.py -i 3.mapped -o 4.stop -g genome.txt -w 8 \\
      --end-type 5 --plus-shift -1 --minus-shift 1 --pattern "*.rmdup.bam"

命令行参数:
  -i, --input-dir     BAM文件输入目录 (默认: 当前目录)
  -o, --output-dir    输出目录 (默认: 当前目录)
  -g, --genome-file   Genome chrNameLength文件路径 (必需)
  -w, --workers       并行线程数 (默认: 1，串行处理)
  -p, --pattern       BAM文件通配符匹配模式 (默认: *.bam)
  --pattern-begin     文件名起始字符筛选
  --pattern-mid       文件名中间部分字符筛选
  -e, --end-type      统计端点类型: 5 或 3 (默认: 5)
  --plus-shift        正链shift值 (默认: 0)
  --minus-shift       负链shift值 (默认: 0)
  -k, --keep-temp     保留中间文件 (默认: 删除)

================================================================================
注意事项:

  1. 确保bedtools已安装并在PATH中
  2. genome文件必须与BAM文件使用相同的染色体命名规则
  3. 日志文件自动保存在输出目录的logs/文件夹下
  4. 并行处理时使用ThreadPoolExecutor，适合I/O密集型任务

================================================================================
"""

import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
from pathlib import Path
import time
from datetime import datetime, timedelta
import logging
from collections import Counter
import pysam

# ---------- 配置日志 ----------
def setup_logger(log_file):
    """设置日志配置"""
    logger = logging.getLogger('bedtools_processing')
    logger.setLevel(logging.INFO)
    
    # 文件处理器
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # 控制台处理器
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # 自定义格式化器，支持毫秒精度（3位小数）
    class MillisecondFormatter(logging.Formatter):
        def formatTime(self, record, datefmt=None):
            ct = self.converter(record.created)
            if datefmt:
                s = time.strftime(datefmt, ct)
            else:
                t = time.strftime('%Y-%m-%d %H:%M:%S', ct)
                s = '%s.%03d' % (t, record.msecs)
            return s
    
    formatter = MillisecondFormatter('[%(asctime)s] %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

# ---------- 工具函数 ----------
def get_duration(start_time, end_time):
    """计算时间差并格式化为 HH:MM:SS.mmm"""
    duration = end_time - start_time
    hours = int(duration // 3600)
    minutes = int((duration % 3600) // 60)
    seconds = duration % 60  # 保留小数部分
    return f"{hours:02d}:{minutes:02d}:{seconds:06.3f}"

def run_command(cmd, check=True):
    """运行shell命令"""
    try:
        result = subprocess.run(cmd, shell=True if isinstance(cmd, str) else False, 
                              capture_output=True, text=True, check=check)
        return result
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"命令执行失败: {' '.join(cmd) if isinstance(cmd, list) else cmd}\n"
                          f"错误信息: {e.stderr}")


def count_endpoints_direct(bam_file, genome_file, output_file, end_type,
                           plus_shift, minus_shift):
    """单次读取 BAM，直接累计链特异性端点，避免多轮 BED 中间文件。"""
    chrom_lengths = {}
    with open(genome_file) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                chrom_lengths[fields[0]] = int(fields[1])

    counts = Counter()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_end is None:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            chrom_len = chrom_lengths.get(chrom)
            if chrom_len is None:
                continue
            strand = '-' if read.is_reverse else '+'
            shift = minus_shift if strand == '-' else plus_shift
            shifted_start = max(0, min(chrom_len, read.reference_start + shift))
            shifted_end = max(0, min(chrom_len, read.reference_end + shift))
            if shifted_start >= shifted_end:
                continue
            if end_type == '5':
                pos = shifted_end - 1 if strand == '-' else shifted_start
            else:
                pos = shifted_start if strand == '-' else shifted_end - 1
            counts[(chrom, pos, strand)] += 1

    with open(output_file, "w") as out:
        for (chrom, pos, strand), count in sorted(
            counts.items(), key=lambda item: (item[0][0], item[0][1], item[0][2] == '-')
        ):
            out.write(f"{chrom}\t{pos}\t{pos + 1}\t{count}\t{strand}\n")
    return len(counts)

# ---------- 核心处理函数 ----------
def process_single_bam(bam_file, genome_file, output_dir, logger, file_index, total_files, 
                       end_type='5', plus_shift=-1, minus_shift=1, keep_temp=False):
    """
    处理单个BAM文件
    
    参数:
      bam_file: BAM文件路径
      genome_file: genome文件路径
      output_dir: 输出目录路径
      logger: 日志对象
      file_index: 当前文件索引
      total_files: 总文件数
      end_type: '5' 或 '3'，指定统计reads的5'端还是3'端
      plus_shift: 正链shift值（默认-1）
      minus_shift: 负链shift值（默认1）
      keep_temp: 是否保留中间文件（默认False）
    """
    try:
        bam_path = Path(bam_file)
        prefix = bam_path.stem  # 获取不带扩展名的文件名
        output_path = Path(output_dir)
        
        # 记录文件处理开始时间
        file_start_time = time.time()
        
        logger.info("")
        logger.info(f"========== 处理文件 [{file_index}/{total_files}]: {bam_path.name} ==========")
        logger.info(f"输出目录: {output_dir}")
        logger.info(f"参数设置: 统计{end_type}'端, 正链shift={plus_shift}, 负链shift={minus_shift}")
        logger.info(f"中间文件: {'保留' if keep_temp else '删除'}")

        combined_sorted = str(output_path / f"{prefix}_shifted_combined_sorted.bedgraph")
        if not keep_temp:
            logger.info("使用单次BAM遍历直接累计端点（无中间BED/bedGraph）...")
            site_count = count_endpoints_direct(
                bam_file, genome_file, combined_sorted, end_type,
                plus_shift, minus_shift
            )
            file_end_time = time.time()
            logger.info(f"端点位点数: {site_count:,}")
            logger.info(f"最终输出: {combined_sorted}")
            logger.info(f"文件处理总耗时: {get_duration(file_start_time, file_end_time)}")
            return True, prefix
        
        # 记录中间文件路径
        temp_files = []
        
        # Step 1: BAM to BED
        step_start = time.time()
        logger.info("开始: BAM转BED...")
        bed_file = str(output_path / f"{prefix}.bed")
        temp_files.append(bed_file)
        cmd = f"bedtools bamtobed -i {bam_file} > {bed_file}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"完成: BAM转BED (耗时: {get_duration(step_start, step_end)})")
        
        # Step 2: Shift
        step_start = time.time()
        logger.info(f"开始: Shift操作 (正链:{plus_shift}, 负链:{minus_shift})...")
        shifted_bed = str(output_path / f"{prefix}_shifted.bed")
        temp_files.append(shifted_bed)
        cmd = f"bedtools shift -m {minus_shift} -p {plus_shift} -i {bed_file} -g {genome_file} > {shifted_bed}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"完成: Shift操作 (耗时: {get_duration(step_start, step_end)})")
        
        # 根据end_type设置genomecov参数
        end_param = f"-{end_type}"  # -5 或 -3
        
        # Step 3: Plus strand coverage
        step_start = time.time()
        logger.info(f"开始: 计算plus链覆盖度 (统计{end_type}'端)...")
        plus_bedgraph = str(output_path / f"{prefix}_shifted_plus.bedgraph")
        temp_files.append(plus_bedgraph)
        cmd = f"bedtools genomecov -bg -strand + {end_param} -i {shifted_bed} -g {genome_file} > {plus_bedgraph}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"完成: Plus链覆盖度 (耗时: {get_duration(step_start, step_end)})")
        
        # Step 4: Minus strand coverage
        step_start = time.time()
        logger.info(f"开始: 计算minus链覆盖度 (统计{end_type}'端)...")
        minus_bedgraph = str(output_path / f"{prefix}_shifted_minus.bedgraph")
        temp_files.append(minus_bedgraph)
        cmd = f"bedtools genomecov -bg -strand - {end_param} -i {shifted_bed} -g {genome_file} > {minus_bedgraph}"
        run_command(cmd)
        step_end = time.time()
        logger.info(f"完成: Minus链覆盖度 (耗时: {get_duration(step_start, step_end)})")
        
        # Step 5: Merge plus and minus strand with strand information
        step_start = time.time()
        logger.info("开始: 合并正负链数据并添加strand信息...")
        
        combined_bedgraph = str(output_path / f"{prefix}_shifted_combined.bedgraph")
        temp_files.append(combined_bedgraph)
        combined_sorted = str(output_path / f"{prefix}_shifted_combined_sorted.bedgraph")
        
        # 读取并处理plus链数据
        with open(plus_bedgraph, 'r') as f_plus, \
             open(combined_bedgraph, 'w') as f_out:
            for line in f_plus:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    f_out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t+\n")
        
        # 追加minus链数据
        with open(minus_bedgraph, 'r') as f_minus, \
             open(combined_bedgraph, 'a') as f_out:
            for line in f_minus:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    f_out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t-\n")
        
        # 排序
        cmd = f"LC_COLLATE=C sort -k1,1 -k2,2n {combined_bedgraph} > {combined_sorted}"
        run_command(cmd)
        
        step_end = time.time()
        logger.info(f"完成: 合并正负链数据 (耗时: {get_duration(step_start, step_end)})")
        
        # 删除中间文件（如果不保留）
        if not keep_temp:
            logger.info("清理: 删除中间文件...")
            for temp_file in temp_files:
                try:
                    if Path(temp_file).exists():
                        Path(temp_file).unlink()
                        logger.debug(f"  已删除: {temp_file}")
                except Exception as e:
                    logger.warning(f"  删除文件 {temp_file} 失败: {e}")
            logger.info("清理: 中间文件已删除")
        else:
            logger.info("保留: 中间文件已保留")
        
        # 计算文件处理总时间
        file_end_time = time.time()
        file_total_duration = get_duration(file_start_time, file_end_time)
        
        logger.info(f"文件 {bam_path.name} 处理完成")
        logger.info(f"最终输出: {combined_sorted}")
        logger.info(f"文件处理总耗时: {file_total_duration}")
        
        return True, prefix
        
    except Exception as e:
        # 即使出错也记录处理时间
        file_end_time = time.time()
        file_total_duration = get_duration(file_start_time, file_end_time)
        logger.error(f"[ERROR] 处理 {bam_file} 时出错: {e}")
        logger.info(f"处理时间 (失败): {file_total_duration}")
        return False, None

def process_bam_parallel(input_dir, output_dir, genome_file, workers=1, pattern='*.bam',
                         pattern_begin=None, pattern_mid=None,
                         end_type='5', plus_shift=-1, minus_shift=1, keep_temp=False):
    """
    并行处理BAM文件
    
    参数:
      input_dir: 输入目录
      output_dir: 输出目录
      genome_file: genome文件路径 (chrNameLength.txt)
      workers: 并行进程数
      pattern: 文件匹配模式 (当pattern_begin和pattern_mid都为None时使用)
      pattern_begin: 文件名起始字符，用于筛选特定样本
      pattern_mid: 文件名中间部分字符，用于筛选特定样本
      end_type: '5' 或 '3'，指定统计reads的5'端还是3'端
      plus_shift: 正链shift值
      minus_shift: 负链shift值
      keep_temp: 是否保留中间文件
    """
    # 确保输出目录存在
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 设置日志 - 输出到输出目录的logs文件夹
    logs_dir = Path(output_dir) / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = str(logs_dir / f"bedtools_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logger = setup_logger(log_file)
    
    # 设置环境变量
    os.environ['LC_COLLATE'] = 'C'
    
    # 记录开始时间
    script_start = time.time()
    logger.info("========== 脚本开始执行 ==========")
    logger.info(f"输入目录: {input_dir}")
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"Genome文件: {genome_file}")
    logger.info(f"参数配置: 统计{end_type}'端, 正链shift={plus_shift}, 负链shift={minus_shift}")
    logger.info(f"中间文件处理: {'保留' if keep_temp else '删除'}")
    
    # 查找BAM文件 - 使用灵活的匹配模式
    input_path = Path(input_dir)
    
    # 优先使用pattern_begin和pattern_mid进行筛选
    if pattern_begin or pattern_mid:
        if pattern_begin:
            if pattern_mid:
                # 同时使用起始和中间匹配
                all_bams = list(input_path.glob(f'{pattern_begin}*.bam'))
                bam_files = [f for f in all_bams if pattern_mid in f.name]
                logger.info(f"使用匹配模式: 起始='{pattern_begin}', 中间='{pattern_mid}'")
            else:
                # 只使用起始匹配
                bam_files = list(input_path.glob(f'{pattern_begin}*.bam'))
                logger.info(f"使用匹配模式: 起始='{pattern_begin}'")
        else:
            # 只使用中间匹配
            all_bams = list(input_path.glob('*.bam'))
            bam_files = [f for f in all_bams if pattern_mid in f.name]
            logger.info(f"使用匹配模式: 中间='{pattern_mid}'")
    else:
        # 使用传统的pattern模式
        bam_files = list(input_path.glob(pattern))
        logger.info(f"使用匹配模式: 通配符='{pattern}'")
    
    if not bam_files:
        if pattern_begin or pattern_mid:
            logger.error(f"在 {input_dir} 下未找到符合条件的BAM文件")
            logger.error(f"匹配条件: 起始='{pattern_begin}', 中间='{pattern_mid}'")
        else:
            logger.error(f"在 {input_dir} 下未找到符合模式 {pattern} 的BAM文件")
        return
    
    bam_count = len(bam_files)
    logger.info(f"找到 {bam_count} 个BAM文件待处理")
    
    # 显示匹配到的文件列表（最多显示10个，超过则省略）
    if bam_count > 0:
        logger.info("匹配到的文件:")
        for i, bam_file in enumerate(bam_files[:10], 1):
            logger.info(f"  {i}. {bam_file.name}")
        if bam_count > 10:
            logger.info(f"  ... 还有 {bam_count - 10} 个文件未显示")
    
    if workers == 1:
        # 串行处理
        for idx, bam_file in enumerate(bam_files, 1):
            process_single_bam(str(bam_file), genome_file, output_dir, logger, idx, bam_count,
                             end_type, plus_shift, minus_shift, keep_temp)
    else:
        # 并行处理
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = []
            for idx, bam_file in enumerate(bam_files, 1):
                futures.append(
                    executor.submit(
                        process_single_bam,
                        str(bam_file),
                        genome_file,
                        output_dir,
                        logger,
                        idx,
                        bam_count,
                        end_type,
                        plus_shift,
                        minus_shift,
                        keep_temp
                    )
                )
            
            # 等待所有任务完成
            success_count = 0
            fail_count = 0
            for future in as_completed(futures):
                try:
                    success, prefix = future.result()
                    if success:
                        success_count += 1
                    else:
                        fail_count += 1
                except Exception as e:
                    logger.error(f"[EXCEPTION] 处理样本时出错: {e}")
                    fail_count += 1
            
            logger.info(f"处理完成 - 成功: {success_count}, 失败: {fail_count}")
    
    # 记录结束时间
    script_end = time.time()
    total_duration = get_duration(script_start, script_end)
    logger.info("")
    logger.info("========== 脚本执行完成 ==========")
    logger.info(f"总耗时: {total_duration}")
    logger.info(f"日志文件: {log_file}")

# ---------- 主程序入口 ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="并行处理BAM文件进行位点鉴定",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法 - 处理当前目录下所有BAM文件
  python call_stop.py
  
  # 指定输入目录、输出目录和genome文件
  python call_stop.py -i /path/to/bam/files -o /path/to/output -g /path/to/genome.txt
  
  # 使用文件名匹配模式并指定输出目录
  python call_stop.py -o ./results --pattern-begin "sample" --pattern-mid "treated"
  python call_stop.py -o ./output --pattern-begin "ctrl" 
  python call_stop.py -o ./processed --pattern-mid "rep1"
  
  # 传统通配符模式（当pattern-begin和pattern-mid都未指定时使用）
  python call_stop.py -o ./results --pattern "*.sorted.bam"
  
  # 并行处理并指定输出目录
  python call_stop.py -w 4 -o ./parallel_output --pattern-begin "sample"
  
  # 自定义参数并指定输出目录
  python call_stop.py -o ./custom_output --end-type 3 --plus-shift 0 --minus-shift 0 --keep-temp
        """)
    parser.add_argument('--input-dir', '-i', default='.', help='BAM文件目录（默认：当前目录）')
    parser.add_argument('--output-dir', '-o', default='.', help='输出文件目录（默认：当前目录）')
    parser.add_argument('--genome-file', '-g', 
                       default='/home/pf/14T/index/human/hg38/human_hg38_genome_star/chrNameLength.txt',
                       help='Genome chrNameLength文件路径')
    parser.add_argument('--workers', '-w', type=int, default=1, 
                       help='并行处理的进程数（默认：1，串行处理）')
    parser.add_argument('--pattern', '-p', default='*.bam',
                       help='BAM文件匹配模式（默认：*.bam，当--pattern-begin和--pattern-mid都未指定时使用）')
    parser.add_argument('--pattern-begin', type=str,
                       help='文件名起始字符，用于筛选特定样本')
    parser.add_argument('--pattern-mid', type=str,
                       help='文件名中间部分字符，用于筛选特定样本')
    parser.add_argument('--end-type', '-e', choices=['5', '3'], default='5',
                       help="统计reads的端点类型：'5'表示5'端，'3'表示3'端（默认：5'端）")
    parser.add_argument('--plus-shift', type=int, default=0,
                       help='正链shift值（默认：0）')
    parser.add_argument('--minus-shift', type=int, default=0,
                       help='负链shift值（默认：0）')
    parser.add_argument('--keep-temp', '-k', action='store_true',
                       help='保留中间文件（默认：删除中间文件，只保留最终结果）')
    
    args = parser.parse_args()
    
    # 检查genome文件是否存在
    if not Path(args.genome_file).exists():
        print(f"[ERROR] Genome文件不存在: {args.genome_file}")
        exit(1)
    
    # 执行处理
    process_bam_parallel(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        genome_file=args.genome_file,
        workers=args.workers,
        pattern=args.pattern,
        pattern_begin=args.pattern_begin,
        pattern_mid=args.pattern_mid,
        end_type=args.end_type,
        plus_shift=args.plus_shift,
        minus_shift=args.minus_shift,
        keep_temp=args.keep_temp
    )
