class LogDataAnalyzer:
    """日志数据分析器，用于提取指定时刻的数据"""
    
    def __init__(self, config=None):
        """
        初始化分析器
        Args:
            config: 配置字典，如果不提供则使用默认配置
        """
        # 默认配置
        self.default_config = {
            'log_file_path': 'log',          # 日志文件路径
            'target_time': 0.6,              # 目标时刻
            'data_prefix': 'pc0[0]:',        # 目标数据前缀
            'time_prefix': 'ti:',            # 时间前缀
            'float_tolerance': 1e-6,         # 浮点数比较容差
            'encoding': 'utf-8',             # 文件编码
        }
        
        # 合并用户配置
        self.config = self.default_config.copy()
        if config:
            self.config.update(config)
        
        # 分析结果
        self.results = {}
    
    def collect_data_at_time(self):
        """
        收集指定时刻的数据并计算统计信息
        Returns:
            dict: 包含统计信息的结果字典
        """
        # 存储所有事件的数据：每个事件是一个列表，包含该事件在目标时刻的所有值
        all_events_data = []
        current_event_data = []
        current_ti = None
        last_ti = float('inf')
        event_count = 0
        
        try:
            with open(self.config['log_file_path'], 'r', 
                     encoding=self.config['encoding']) as file:
                for line_num, line in enumerate(file, 1):
                    line = line.strip()
                    
                    # 检测时间行
                    if line.startswith(self.config['time_prefix']):
                        try:
                            ti_value = float(line.split(':')[1].strip())
                            
                            # 检查是否是新事件的开始（ti变小）
                            if ti_value < last_ti:
                                event_count += 1
                                # 处理上一个事件的数据
                                if current_event_data:
                                    all_events_data.append(current_event_data)
                                    current_event_data = []
                            
                            current_ti = ti_value
                            last_ti = ti_value
                            
                        except ValueError as e:
                            print(f"警告: 第{line_num}行时间格式错误: {line}")
                            continue
                    
                    # 检测目标数据行
                    elif line.startswith(self.config['data_prefix']):
                        try:
                            value = float(line.split(':')[1].strip())
                            
                            if current_ti is not None:
                                # 检查是否为目标时刻
                                if abs(current_ti - self.config['target_time']) < self.config['float_tolerance']:
                                    current_event_data.append(value)
                                    
                        except ValueError as e:
                            print(f"警告: 第{line_num}行数据格式错误: {line}")
                            continue
            
            # 处理最后一个事件的数据
            if current_event_data:
                all_events_data.append(current_event_data)
        
        except FileNotFoundError:
            print(f"错误: 找不到文件 '{self.config['log_file_path']}'")
            return self._create_empty_result()
        except Exception as e:
            print(f"错误: {e}")
            return self._create_empty_result()
        
        # 计算统计信息
        return self._calculate_statistics(all_events_data, event_count)
    
    def _calculate_statistics(self, all_events_data, event_count):
        """计算统计信息"""
        # 将所有数据展平为一个列表
        all_values = []
        for event_data in all_events_data:
            all_values.extend(event_data)
        
        if not all_values:
            return self._create_empty_result()
        
        # 基础统计
        avg_value = sum(all_values) / len(all_values)
        min_value = min(all_values)
        max_value = max(all_values)
        
        # 计算标准差
        if len(all_values) > 1:
            variance = sum((x - avg_value) ** 2 for x in all_values) / len(all_values)
            std_dev = variance ** 0.5
        else:
            std_dev = 0
        
        # 计算有数据的事件数
        events_with_data = len(all_events_data)
        
        # 构建结果字典
        result = {
            'target_time': self.config['target_time'],
            'data_type': self.config['data_prefix'].rstrip(':'),
            'total_events': event_count,
            'events_with_data': events_with_data,
            'all_events_data': all_events_data,  # 每个事件的数据列表
            'all_values': all_values,           # 所有值的扁平列表
            'average': avg_value,
            'min': min_value,
            'max': max_value,
            'std_dev': std_dev,
            'range': max_value - min_value,
            'success_rate': events_with_data / event_count if event_count > 0 else 0,
            'total_data_points': len(all_values),  # 总数据点数
        }
        
        self.results = result
        return result
    
    def _create_empty_result(self):
        """创建空结果字典"""
        return {
            'target_time': self.config['target_time'],
            'data_type': self.config['data_prefix'].rstrip(':'),
            'total_events': 0,
            'events_with_data': 0,
            'all_events_data': [],
            'all_values': [],
            'average': None,
            'min': None,
            'max': None,
            'std_dev': None,
            'range': None,
            'success_rate': 0,
            'total_data_points': 0,
        }
    
    def print_results(self):
        """打印分析结果（统计信息在最后）"""
        if not self.results:
            print("没有可用的分析结果")
            return
        
        print("=" * 60)
        print(f"数据收集报告 - {self.results['data_type']} @ ti={self.results['target_time']}")
        print("=" * 60)
        
        # 1. 首先打印所有收集到的数据
        print(f"\n收集到的数据值（每个事件一个组）:")
        print("-" * 60)
        
        if self.results['all_events_data']:
            for i, event_data in enumerate(self.results['all_events_data'], 1):
                if event_data:  # 只显示有数据的事件
                    print(f"事件 {i}:")
                    for j, value in enumerate(event_data, 1):
                        print(f"    数据点 {j}: {value:.8f}")
                    # 如果事件有多个数据点，显示该事件的统计
                    if len(event_data) > 1:
                        event_avg = sum(event_data) / len(event_data)
                        print(f"    本事件平均: {event_avg:.8f} (共{len(event_data)}个数据点)")
                    print()
        else:
            print("    未收集到任何数据")
        
        # 2. 然后打印统计信息
        print("-" * 60)
        print(f"统计摘要:")
        print("-" * 60)
        print(f"总事件数: {self.results['total_events']}")
        print(f"包含数据的事件数: {self.results['events_with_data']}")
        print(f"数据收集成功率: {self.results['success_rate']:.2%}")
        print(f"总数据点数: {self.results['total_data_points']}")
        print("-" * 60)
        print(f"总体平均值: {self.results['average']:.8f}")
        print(f"最小值: {self.results['min']:.8f}")
        print(f"最大值: {self.results['max']:.8f}")
        print(f"标准差: {self.results['std_dev']:.8f}")
        print(f"数据范围: {self.results['range']:.8f}")
        print("=" * 60)


# 配置和使用示例
if __name__ == "__main__":
    # ========== 参数配置区域 ==========
    # 在这里修改配置参数
    
    # 单个分析配置
    SINGLE_ANALYSIS_CONFIG = {
        'log_file_path': 'log',      # 日志文件路径
        'target_time': 4.2,          # 目标时刻
        'data_prefix': 'np:',    # 要收集的数据前缀（可改为'mass:'或'totprob:'等）
        'time_prefix': 'ti:',        # 时间标记前缀
        'float_tolerance': 1e-6,     # 浮点数比较容差
        'encoding': 'utf-8',         # 文件编码
    }
    
    # 批量分析多个时间点
    MULTIPLE_TIMES_CONFIG = {
        'log_file_path': 'log',
        'data_prefix': 'energy:',
        'target_times': [3.0,5.1,7.2],  # 要分析的多个时间点
    }
    
    # 批量分析不同类型的数据
    MULTIPLE_DATA_TYPES_CONFIG = {
        'log_file_path': 'log',
        'target_time': 7.2,
        'data_types': ['pc0[0]:',  'prob23:','pc01[0]:','pc01new[0]:','fg:'],  # 要分析的多个数据类型
    }
    # ========== 配置结束 ==========
    
    # 方法1：分析单个时间点和数据类型
    print("单个时间点分析:")
    analyzer = LogDataAnalyzer(SINGLE_ANALYSIS_CONFIG)
    result = analyzer.collect_data_at_time()
    analyzer.print_results()
    
    print("\n\n")
    '''
    # 方法2：批量分析多个时间点（可选）
    print("多个时间点批量分析:")
    print("=" * 60)
    
    for time_point in MULTIPLE_TIMES_CONFIG['target_times']:
        config = {
            'log_file_path': MULTIPLE_TIMES_CONFIG['log_file_path'],
            'target_time': time_point,
            'data_prefix': MULTIPLE_TIMES_CONFIG['data_prefix'],
        }
        
        analyzer = LogDataAnalyzer(config)
        result = analyzer.collect_data_at_time()
        
        if result['total_data_points'] > 0:
            print(f"时刻 {time_point:.1f}: "
                  f"平均值={result['average']:.8f}, "
                  f"数据点={result['total_data_points']}, "
                  f"事件数={result['events_with_data']}/{result['total_events']}")
        else:
            print(f"时刻 {time_point:.1f}: 未找到数据")
    
    print("=" * 60)
    '''
    # 方法3：分析不同类型的数据（可选）
    print("\n\n分析不同类型的数据:")
    print("=" * 60)
    
    for data_type in MULTIPLE_DATA_TYPES_CONFIG['data_types']:
        config = {
            'log_file_path': MULTIPLE_DATA_TYPES_CONFIG['log_file_path'],
            'target_time': MULTIPLE_DATA_TYPES_CONFIG['target_time'],
            'data_prefix': data_type,
        }
        
        analyzer = LogDataAnalyzer(config)
        result = analyzer.collect_data_at_time()
        
        if result['total_data_points'] > 0:
            print(f"\n{data_type.rstrip(':')} (ti={MULTIPLE_DATA_TYPES_CONFIG['target_time']}):")
            print(f"    平均值={result['average']:.8f}, "
                  f"数据点={result['total_data_points']}, "
                  f"事件数={result['events_with_data']}/{result['total_events']}")
        else:
            print(f"\n{data_type.rstrip(':')} (ti={MULTIPLE_DATA_TYPES_CONFIG['target_time']}): 未找到数据")
    
    print("=" * 60)
    
