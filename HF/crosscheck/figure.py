import numpy as np
import matplotlib.pyplot as plt
import math

# 定义计算φ的函数
def psi_std(px, py):
    if np.isnan(px):
        px = 0
    if np.isnan(py):
        py = 0

    if px == 0 and py == 0:
        psi = 0
    elif px == 0 and py > 0:
        psi = math.pi / 2
    elif px == 0 and py < 0:
        psi = -math.pi / 2
    elif px > 0:
        psi = math.atan(py / px)
    elif px < 0 and py >= 0:
        psi = math.pi + math.atan(py / px)
    elif px < 0 and py < 0:
        psi = -math.pi + math.atan(py / px)
    return psi

# 读取第一个文件
def read_file1(filename):
    data = np.loadtxt(filename)
    event_numbers = data[:, 0]
    pdg_numbers = data[:, 1]
    px = data[:, 2]
    py = data[:, 3]
    pz = data[:, 4]
    E = data[:, 5]
    return event_numbers, pdg_numbers, px, py, pz, E

# 读取第二个文件
def read_file2(filename):
    data = np.loadtxt(filename)
    event_numbers = data[:, 0]
    pdg_numbers = data[:, 1]
    px = data[:, 2]
    py = data[:, 3]
    pz = data[:, 4]
    E = data[:, 5]
    return event_numbers, pdg_numbers, px, py, pz, E

# 计算物理量
def calculate_quantities(px, py, pz, E):
    # 计算pT
    pt = np.sqrt(px**2 + py**2)
    
    # 计算η
    p = np.sqrt(px**2 + py**2 + pz**2)
    # 避免除零错误
    mask = (p != np.abs(pz))
    eta = np.zeros_like(pz)
    eta[mask] = 0.5 * np.log((p[mask] + pz[mask]) / (p[mask] - pz[mask]))
    
    # 计算φ
    phi = np.array([psi_std(px_i, py_i) for px_i, py_i in zip(px, py)])
    
    return pt, eta, phi, E

# 绘制分布图
def plot_distributions(data1, data2, num_events1, num_events2, bins=50):
    pt1, eta1, phi1, E1 = data1
    pt2, eta2, phi2, E2 = data2
    
    # 创建图形
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # dN/dE分布
    hist_E1, bin_edges_E = np.histogram(E1, bins=bins)
    hist_E2, _ = np.histogram(E2, bins=bin_edges_E)
    
    bin_centers_E = (bin_edges_E[:-1] + bin_edges_E[1:]) / 2
    bin_widths_E = np.diff(bin_edges_E)
    
    axes[0, 0].plot(bin_centers_E, hist_E1 / (num_events1 * bin_widths_E), 
                   label='PYTHIA+Frag', linewidth=2)
    axes[0, 0].plot(bin_centers_E, hist_E2 / (num_events2 * bin_widths_E), 
                   label='PYTHIA', linewidth=2)
    axes[0, 0].set_xlabel('E (GeV)')
    axes[0, 0].set_ylabel('dN/dE')
    axes[0, 0].legend()
    axes[0, 0].grid(True)
    
    # dN/dφ分布
    hist_phi1, bin_edges_phi = np.histogram(phi1, bins=bins, range=(-np.pi, np.pi))
    hist_phi2, _ = np.histogram(phi2, bins=bin_edges_phi)
    
    bin_centers_phi = (bin_edges_phi[:-1] + bin_edges_phi[1:]) / 2
    bin_widths_phi = np.diff(bin_edges_phi)
    
    axes[0, 1].plot(bin_centers_phi, hist_phi1 / (num_events1 * bin_widths_phi), 
                   label='PYTHIA+Frag', linewidth=2)
    axes[0, 1].plot(bin_centers_phi, hist_phi2 / (num_events2 * bin_widths_phi), 
                   label='PYTHIA', linewidth=2)
    axes[0, 1].set_xlabel('φ (rad)')
    axes[0, 1].set_ylabel('dN/dφ')
    axes[0, 1].legend()
    axes[0, 1].grid(True)
    
    # dN/dη分布
    # 限制η的范围以避免无穷大值
    eta_range = (-5, 5)
    hist_eta1, bin_edges_eta = np.histogram(eta1, bins=bins, range=eta_range)
    hist_eta2, _ = np.histogram(eta2, bins=bin_edges_eta)
    
    bin_centers_eta = (bin_edges_eta[:-1] + bin_edges_eta[1:]) / 2
    bin_widths_eta = np.diff(bin_edges_eta)
    
    axes[1, 0].plot(bin_centers_eta, hist_eta1 / (num_events1 * bin_widths_eta), 
                   label='PYTHIA+Frag', linewidth=2)
    axes[1, 0].plot(bin_centers_eta, hist_eta2 / (num_events2 * bin_widths_eta), 
                   label='PYTHIA', linewidth=2)
    axes[1, 0].set_xlabel('η')
    axes[1, 0].set_ylabel('dN/dη')
    axes[1, 0].legend()
    axes[1, 0].grid(True)
    
    # dN/dpT分布
    pt_range = (0, 10)  # 假设pT最大为10 GeV
    hist_pt1, bin_edges_pt = np.histogram(pt1, bins=bins, range=pt_range)
    hist_pt2, _ = np.histogram(pt2, bins=bin_edges_pt)
    
    bin_centers_pt = (bin_edges_pt[:-1] + bin_edges_pt[1:]) / 2
    bin_widths_pt = np.diff(bin_edges_pt)
    
    axes[1, 1].plot(bin_centers_pt, hist_pt1 / (num_events1 * bin_widths_pt), 
                   label='PYTHIA+Frag', linewidth=2)
    axes[1, 1].plot(bin_centers_pt, hist_pt2 / (num_events2 * bin_widths_pt), 
                   label='PYTHIA', linewidth=2)
    axes[1, 1].set_xlabel('pT (GeV)')
    axes[1, 1].set_ylabel('dN/dpT')
    axes[1, 1].legend()
    axes[1, 1].grid(True)
    
    plt.tight_layout()
    plt.savefig("show.png")

# 主函数
def main():
    # 读取数据文件
    file1 = '../../result/frag2.dat'
    file2 = './outputdatafile/check.dat'  # 请替换为实际文件名
    
    # 读取第一个文件
    event_numbers1, _, px1, py1, pz1, E1 = read_file1(file1)
    num_events1 = 5000
    
    # 读取第二个文件
    event_numbers2, _, px2, py2, pz2, E2 = read_file2(file2)
    num_events2 = 5000
    
    # 计算物理量
    pt1, eta1, phi1, E1 = calculate_quantities(px1, py1, pz1, E1)
    pt2, eta2, phi2, E2 = calculate_quantities(px2, py2, pz2, E2)
    
    # 绘制分布图
    data1 = (pt1, eta1, phi1, E1)
    data2 = (pt2, eta2, phi2, E2)
    
    plot_distributions(data1, data2, num_events1, num_events2)

if __name__ == '__main__':
    main()
