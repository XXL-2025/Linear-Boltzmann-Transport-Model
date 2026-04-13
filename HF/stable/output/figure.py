import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_dnde_log(data_file, initial_energy=100.0, de_max=100.0):
    """
    绘制dn/dΔE vs ΔE的对数阶梯图
    
    参数:
    data_file: 数据文件名
    initial_energy: 初始能量 (GeV)
    de_max: ΔE的最大值 (GeV)
    """

    try:
        # 读取数据文件
        print(f"Reading data from: {data_file}")

        # 读取数据（假设第五列是能量）
        data = np.loadtxt(data_file)

        # 获取第五列（索引为4）的能量数据
        if data.ndim == 1:
            energies = data  # 如果只有一列
        else:
            energies = data[:, 4]  # 第五列
        
        print(f"Loaded {len(energies)} particles")
        print(f"Energy range: {energies.min():.2f} - {energies.max():.2f} GeV")

        # 计算能量损失 ΔE = 初始能量 - 最终能量
        energy_loss = initial_energy - energies

        # 确保ΔE在合理范围内
        valid_mask = (energy_loss >= 0) & (energy_loss <= de_max)
        energy_loss = energy_loss[valid_mask]

        print(f"Valid data points: {len(energy_loss)}")
        print(f"ΔE range: {energy_loss.min():.2f} - {energy_loss.max():.2f} GeV")

        # 创建直方图
        bins = 100  # 使用100个bin，每个1 GeV
        counts, bin_edges = np.histogram(energy_loss, bins=bins, range=(0, de_max))

        # 计算bin宽度
        bin_width = bin_edges[1] - bin_edges[0]

        # 计算dn/dΔE = counts / (bin_width * N_total)
        n_total = len(energy_loss)
        dnde = counts / (bin_width * n_total)

        # 计算bin中心
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # 创建图形
        plt.figure(figsize=(10, 7))

        # 绘制dn/dΔE vs ΔE阶梯图（对数纵坐标）
        plt.step(bin_centers, dnde, where='mid', linewidth=2, color='blue')

        # 设置对数纵坐标
       # plt.yscale('log')

        # 设置轴标签
        plt.xlabel('ΔE (GeV)', fontsize=14)
        plt.ylabel('dn/dΔE (1/GeV)', fontsize=14)

        # 设置网格
        plt.grid(True, alpha=0.3, which='both')

        # 设置x轴范围
        plt.xlim(0, de_max)

        # 添加统计信息
        mean_de = np.mean(energy_loss)
        median_de = np.median(energy_loss)



        # 调整布局
        plt.tight_layout()
        # 保存图像
        output_file = 'dnde_log.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved as: {output_file}")

        # 显示图像
        plt.show()

        # 打印统计信息
        print("\n===== Statistics =====")
        print(f"Initial energy: {initial_energy} GeV")
        print(f"Number of particles: {n_total}")
        print(f"Mean energy loss: {mean_de:.4f} GeV")
        print(f"Median energy loss: {median_de:.4f} GeV")
        print(f"Std of energy loss: {np.std(energy_loss):.4f} GeV")
        print(f"Mean remaining energy: {np.mean(energies[valid_mask]):.4f} GeV")

        # 计算不同能量损失区间的比例
        de_ranges = [(0, 10), (10, 30), (30, 70), (70, 100)]
        print("\nEnergy loss distribution:")
        for de_min, de_max_range in de_ranges:
            ratio = np.sum((energy_loss >= de_min) & (energy_loss < de_max_range)) / n_total * 100
            print(f"  ΔE in {de_min}-{de_max_range} GeV: {ratio:.2f}%")

        return energy_loss, dnde, bin_centers

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None

# 如果是直接运行脚本
if __name__ == "__main__":
    # 默认参数
    data_file = "1.dat"  # 默认文件名
    initial_energy = 100.0  # GeV
    de_max = 100.0  # GeV
    # 如果提供了命令行参数
    if len(sys.argv) > 1:
        data_file = sys.argv[1]

    # 绘制图像
    plot_dnde_log(data_file, initial_energy, de_max)
