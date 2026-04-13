import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def plot_dnde_log_compare(data_file_hf, data_file_merge, data_file_stable, initial_energy=100.0, de_max=100.0):
    """
    绘制三条dn/dΔE vs ΔE的对数阶梯图进行对比
    
    参数:
    data_file_hf: LBT_HF数据文件名
    data_file_merge: LBT_MERGE数据文件名
    data_file_stable: LBT_HF_STABLE数据文件名
    initial_energy: 初始能量 (GeV) - 用于LBT_HF和LBT_HF_STABLE格式
    de_max: ΔE的最大值 (GeV)
    """
    
    try:
        # 设置图形
        plt.figure(figsize=(10, 7))
        
        # 定义颜色、标签和线型
        colors = ['blue', 'red', 'green']
        labels = ['LBT_HF', 'LBT_MERGE', 'LBT_HF_STABLE']
        line_styles = ['-', '-', '-']
        
        all_energy_loss_data = []
        all_dnde_data = []
        all_bin_centers = []
        
        # 获取要处理的文件列表
        data_files = [data_file_hf, data_file_merge, data_file_stable]
        
        # 处理每个文件
        for i, data_file in enumerate(data_files):
            print(f"\n{'='*50}")
            print(f"Processing {labels[i]}: {data_file}")
            
            # 检查文件是否存在
            if not os.path.exists(data_file):
                print(f"Warning: File {data_file} not found. Skipping...")
                all_energy_loss_data.append(np.array([]))
                all_dnde_data.append(np.array([]))
                all_bin_centers.append(np.array([]))
                continue
            
            # 读取数据
            data = np.loadtxt(data_file)
            
            if i == 0 or i == 2:  # LBT_HF和LBT_HF_STABLE格式
                if data.ndim == 1:
                    energies = data
                else:
                    energies = data[:, 4]  # 第五列是能量
                
                # 计算能量损失 ΔE = 初始能量 - 最终能量
                energy_loss = initial_energy - energies
                
            else:  # LBT_MERGE格式 (i == 1)
                if data.ndim == 1:
                    print("Error: LBT_MERGE data should have at least 3 columns")
                    all_energy_loss_data.append(np.array([]))
                    all_dnde_data.append(np.array([]))
                    all_bin_centers.append(np.array([]))
                    continue
                else:
                    initial_energies = data[:, 1]  # 第二列是初始能量
                    final_energies = data[:, 2]    # 第三列是最终能量
                
                # 计算能量损失 ΔE = 初始能量 - 最终能量
                energy_loss = initial_energies - final_energies
            
            # 确保ΔE在合理范围内
            valid_mask = (energy_loss >= 0) & (energy_loss <= de_max)
            energy_loss = energy_loss[valid_mask]
            
            print(f"Loaded {len(energy_loss)} valid particles")
            if len(energy_loss) > 0:
                print(f"ΔE range: {energy_loss.min():.2f} - {energy_loss.max():.2f} GeV")
            else:
                print("No valid particles found")
            
            # 创建直方图
            bins = 100
            counts, bin_edges = np.histogram(energy_loss, bins=bins, range=(0, de_max))
            
            # 计算bin宽度
            bin_width = bin_edges[1] - bin_edges[0]
            
            # 计算dn/dΔE = counts / (bin_width * N_total)
            n_total = len(energy_loss)
            if n_total > 0:
                dnde = counts / (bin_width * n_total)
            else:
                dnde = np.zeros_like(counts, dtype=float)
            
            # 计算bin中心
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            all_energy_loss_data.append(energy_loss)
            all_dnde_data.append(dnde)
            all_bin_centers.append(bin_centers)
        
        # 绘制所有曲线
        for i in range(len(data_files)):
            if len(all_energy_loss_data[i]) > 0:
                # 绘制阶梯图
                plt.step(all_bin_centers[i], all_dnde_data[i], where='mid', 
                        linewidth=2.5, color=colors[i], linestyle=line_styles[i], 
                        label=labels[i], alpha=0.9)
        
        # 设置图形属性
        #plt.yscale('log')
        plt.xlabel('ΔE (GeV)', fontsize=14)
        plt.ylabel('dn/dΔE (1/GeV)', fontsize=14)
        plt.title('Comparison of Energy Loss Distributions', fontsize=16)
        plt.grid(True, alpha=0.3, which='both')
        plt.xlim(0, de_max)
        plt.ylim(bottom=0)  # 确保y轴从0开始
        plt.legend(fontsize=12, loc='best')
        
        # 调整布局
        plt.tight_layout()
        
        # 保存图像
        output_file = 'dnde_log_comparison_three.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n{'='*50}")
        print(f"Figure saved as: {output_file}")
        
        # 显示图像
        plt.show()
        
        # 打印每个数据集的统计信息
        print(f"\n{'='*50}")
        print("STATISTICS SUMMARY")
        print('='*50)
        
        for i in range(len(data_files)):
            energy_loss = all_energy_loss_data[i]
            if len(energy_loss) > 0:
                print(f"\n{labels[i]}:")
                print(f"  Number of particles: {len(energy_loss)}")
                print(f"  Mean energy loss: {np.mean(energy_loss):.4f} GeV")
                print(f"  Median energy loss: {np.median(energy_loss):.4f} GeV")
                print(f"  Std of energy loss: {np.std(energy_loss):.4f} GeV")
                
                # 计算不同能量损失区间的比例
                de_ranges = [(0, 10), (10, 30), (30, 70), (70, 100)]
                print("  Energy loss distribution:")
                for de_min, de_max_range in de_ranges:
                    ratio = np.sum((energy_loss >= de_min) & (energy_loss < de_max_range)) / len(energy_loss) * 100
                    print(f"    ΔE in {de_min}-{de_max_range} GeV: {ratio:.2f}%")
            else:
                print(f"\n{labels[i]}: No data available")
        
        return all_energy_loss_data, all_dnde_data, all_bin_centers
    
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None

# 单独绘制某个格式的函数（可选）
def plot_dnde_log_single(data_file, label="", initial_energy=100.0, de_max=100.0, format_type="hf"):
    """
    绘制单个数据文件的dn/dΔE vs ΔE的对数阶梯图
    
    参数:
    data_file: 数据文件名
    label: 曲线标签
    initial_energy: 初始能量 (GeV) - 用于hf格式
    de_max: ΔE的最大值 (GeV)
    format_type: 数据格式类型 ("hf" 或 "merge")
    """
    
    try:
        print(f"\n{'='*50}")
        print(f"Processing {label}: {data_file}")
        print('='*50)
        
        # 检查文件是否存在
        if not os.path.exists(data_file):
            print(f"Error: File {data_file} not found.")
            return None, None, None
        
        # 读取数据
        data = np.loadtxt(data_file)
        
        if format_type.lower() == "hf":
            if data.ndim == 1:
                energies = data
            else:
                energies = data[:, 4]  # 第五列是能量
            
            # 计算能量损失
            energy_loss = initial_energy - energies
            
        elif format_type.lower() == "merge":
            if data.ndim == 1:
                print("Error: LBT_MERGE data should have at least 3 columns")
                return None, None, None
            else:
                initial_energies = data[:, 1]  # 第二列是初始能量
                final_energies = data[:, 2]    # 第三列是最终能量
            
            # 计算能量损失
            energy_loss = initial_energies - final_energies
        
        else:
            print(f"Error: Unknown format type {format_type}")
            return None, None, None
        
        # 确保ΔE在合理范围内
        valid_mask = (energy_loss >= 0) & (energy_loss <= de_max)
        energy_loss = energy_loss[valid_mask]
        
        print(f"Loaded {len(energy_loss)} valid particles")
        if len(energy_loss) > 0:
            print(f"ΔE range: {energy_loss.min():.2f} - {energy_loss.max():.2f} GeV")
        
        # 创建直方图
        bins = 100
        counts, bin_edges = np.histogram(energy_loss, bins=bins, range=(0, de_max))
        
        # 计算dn/dΔE
        bin_width = bin_edges[1] - bin_edges[0]
        n_total = len(energy_loss)
        
        if n_total > 0:
            dnde = counts / (bin_width * n_total)
        else:
            dnde = np.zeros_like(counts, dtype=float)
            print("Warning: No valid particles found")
        
        # 计算bin中心
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # 创建图形
        plt.figure(figsize=(10, 7))
        
        # 绘制阶梯图
        plt.step(bin_centers, dnde, where='mid', linewidth=2.5, 
                color='blue', label=label, alpha=0.9)
        
        # 设置图形属性
        #plt.yscale('log')
        plt.xlabel('ΔE (GeV)', fontsize=14)
        plt.ylabel('dn/dΔE (1/GeV)', fontsize=14)
        plt.title(f'Energy Loss Distribution ({label})', fontsize=16)
        plt.grid(True, alpha=0.3, which='both')
        plt.xlim(0, de_max)
        plt.ylim(bottom=0)
        if label:
            plt.legend(fontsize=12)
        
        # 调整布局
        plt.tight_layout()
        
        # 保存图像
        output_file = f'dnde_log_{label.lower().replace(" ", "_")}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved as: {output_file}")
        
        # 显示图像
        plt.show()
        
        # 打印统计信息
        if len(energy_loss) > 0:
            print(f"\n===== {label} Statistics =====")
            print(f"Number of particles: {n_total}")
            print(f"Mean energy loss: {np.mean(energy_loss):.4f} GeV")
            print(f"Median energy loss: {np.median(energy_loss):.4f} GeV")
            print(f"Std of energy loss: {np.std(energy_loss):.4f} GeV")
        
        return energy_loss, dnde, bin_centers
    
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None

# 如果是直接运行脚本
if __name__ == "__main__":
    # 默认参数
    data_file_hf = "save.dat"  # LBT_HF格式文件
    data_file_merge = "/home/hyj/lxx/LBT_MERGE_version2/1.dat"  # LBT_MERGE格式文件
    data_file_stable = "../stable/LBT_HF/output/1.dat"  # LBT_HF_STABLE格式文件
    initial_energy = 100.0  # GeV (用于LBT_HF和LBT_HF_STABLE)
    de_max = 100.0  # GeV
    
    # 如果提供了命令行参数
    if len(sys.argv) >= 4:
        data_file_hf = sys.argv[1]
        data_file_merge = sys.argv[2]
        data_file_stable = sys.argv[3]
    elif len(sys.argv) >= 2:
        print("Usage: python script.py <LBT_HF_file> <LBT_MERGE_file> <LBT_HF_STABLE_file>")
        print("Using default file names")
    
    # 选择要运行的功能：
    
    # 1. 同时绘制三条曲线进行对比（推荐）
    print("\n" + "="*50)
    print("PLOTTING COMPARISON OF ALL THREE DATASETS")
    print("="*50)
    plot_dnde_log_compare(data_file_hf, data_file_merge, data_file_stable, initial_energy, de_max)
    
