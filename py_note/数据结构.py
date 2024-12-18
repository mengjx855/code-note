#### 数据结构 ####
# https://blog.csdn.net/cainiao_python/article/details/111878000


# 堆栈：栈是一个非常重要的数据结构，支持快速后进/先出（LIFO）语义插入和删除。与列表或数组不同，堆栈通常不允许随机访问它们包含的对象。
# 插入和删除操作通常也称为push和pop，Python的栈使用起来要比C,C++简单很多。
# list：
x = []
x.append('eat'), x.append('sleep'), x.append('code')
# ['eat', 'sleep', 'code']
x.pop()
# 'code'
x
# ['eat', 'sleep']
x.append('jump')
# ['eat', 'sleep', 'jump']
   
# deque：双向队列 Python的deque对象以双向链接列表的形式实现,它的操作很像list 同时 ,相比于list实现的队列，deque实现拥有更低的时间和空间复杂度。
# list实现在出队（pop）和插入（insert）时的空间复杂度大约为O(n)，deque在出队（pop）和入队（append）时的时间复杂度是O(1)
# 双向队列其实有很多种妙用，尤其是在访问队列里面的头尾的数据效率非常高，而且还支持很多队列的黑科技，比如逆时针旋转。
from collections import deque
q = deque()
q.append('eat'), q.append('sleep'), q.append('code')
# deque(['eat', 'sleep', 'code'])
q.appendleft('apple')
# deque(['apple', 'eat', 'sleep', 'code'])
q.popleft()
# 'apple'
q
# deque(['eat', 'sleep', 'code'])

# queue 并发里面的队列。当我们处理大量的数据请求的时候，比如我们需要爬大量的网站的网址，有的时候我们会将待处理的请求扔到队列queue里面，用多进程或者多线程进行并发处理。比如典型的生产者消费者的模式中就经常用到queue。
# queue 是一个模块，提供多种类型的队列实现，主要用于多线程编程中。队列是一种先进先出（FIFO, First In First Out）的数据结构，非常适合于任务调度和线程间通信。
# 先进先出：队列中最早的元素最先被移除。线程安全：queue 模块中的队列是为多线程环境设计的，可以在多个线程之间安全地共享。阻塞操作：提供了阻塞队列，当队列为空时，get() 操作会阻塞直到有元素被添加；类似地，当队列满时，put() 操作会阻塞直到有空间可用。
# queue.Queue 还提供了阻塞版本的 get() 和 put() 方法，当队列为空或已满时，这些方法会阻塞调用线程。需要了在研究。
# 几种类型的queue
# Queue：一个标准的先进先出队列，使用锁或其他同步机制来确保线程安全。
# LifoQueue：一个后进先出（LIFO）队列，类似于栈。
# PriorityQueue：一个支持优先级排序的队列，元素根据给定的优先级顺序出队。
# SimpleQueue：一个简单的队列，不提供全部的线程安全特性，但可以用于与不受信任的进程间通信。
import queue
q = queue.Queue() # 创建一个队列实例
q.put('任务1'), q.put('任务2'), q.put('任务3') # 向队列中添加元素
# 从队列中移除并返回元素
while not q.empty():
    task = q.get()  # 获取队列中的元素
    print(task)    # 处理任务
    q.task_done()  # 通知队列该任务已处理完成
# 任务1, 任务2, 任务3




# 字典 Dict：映射(Map)和哈希表(Hash Table)是Python里面的核心数据结构，类似的这样的数据结构几乎在所有的主流的编程语言比如Java ，C++里面也有它的身影
# 普通字典：常见的那种字典的格式
squares = {x:x*x for x in range(6)}

# 顺序字典：在某些情况下，需要保持字典的顺序，这个时候需要用到OrderedDict
# 在 Python 3.7 及以后的版本中，普通字典（dict）在插入元素时会保持插入顺序。这意味着你插入元素的顺序即为元素的存储顺序。在此之前的版本，dict 的顺序是不确定的。
from collections import OrderedDict # Collections 模块提供了一系列的容器数据类型，这些类型是对内置容器数据类型（如 dict, list, set, 和 tuple）的补充和扩展。

# 缺省字典：如果你在存储数据的时候，希望有默认的值在里面，那么你就应该用defaultdict，它在很多场景下面非常巧妙的用处，可以节省很多代码。
from collections import defaultdict
dd_list = defaultdict(list) # 使用list作为默认值类型，也可以是其他类型，例如tuple，dict
dd_list['key1'].append('value1') # 输出：{'key1': ['value1']}
print(dd_list['key2']) # 访问不存在的键，输出：[] 而不是抛出KeyError

# 链式映射：当我们遇到复杂的数据结构的时候，比如需要把多个字典进行合并成一个单一的字典，进行查找搜寻就需要用ChainMap
from collections import ChainMap
d1, d2 = {'one': 1}, {'two': 2}
chain = ChainMap(d1, d2) # 可以随意在chain中访问'one'和'two'




# 数组 array：数组将信息存储在相邻的内存块中，所以被认为是连续的数据结构，很多静态语言里面，都是要强制初始化数组的类型以及长度，如果数组越界会出现寻址到乱的数据
# list可变动态数组
arr = ['one', 'two', 'three']
del arr[1]
arr.append(10)

# tuple不可变数组：跟list一样，但是tuple对象是不可变的，不能修改，删除
arr = ('one', 'two', 'three')

# 类C的array：array模块，可以创建跟C语言很类似的数组，array.array用法跟list很像，唯一的区别在于它只能存储同样地数据类型的数据。它所占的存储空间的大小就是数据的大小。
import array
# 创建 array.array 的实例时，需要指定数组中元素的类型。类型由单字符代码指定，
# b 表示有符号字符，B 表示无符号字符，u 表示无符号短整型，i 表示有符号整型，l 表示有符号长整型，f 表示浮点型，d 表示双精度浮点型
arr_int = array.array('i', [1, 2, 3, 4, 5]) # 创建一个整型数组
arr_float = array.array('d', [1.0, 2.5, 3.5]) # 创建一个浮点型数组

# 字符型数组：使用str对象将文本数据存储为Unicode字符的不可变序列。这意味着 str型字符串数组是不可变的字符数组。
arr = 'abcde'




# 集合：Python也有实现可变且不可变的集合和多集合数据结构。
# set：集合。用来处理去掉重复元素非常有效，多个集合可以进行运算处理，比如取两个集合的交集，联合等等。
x = {'h', 'e', 'l', 'l', 'o'}
'e' in x
x.add('x')
set('Hey') & x
set('Hey') - x

# frozenset 不可变集合 frozenset类实现的不可变set.frozenset对象是静态的，并且仅允许对其元素进行查询操作，而不能进行插入或删除操作。
x = frozenset({'h', 'e', 'l', 'l', 'o'})
x.add('x') # 不可执行

# 多集合Counter：Python标准库中的collection里面的Counter类实现了一种多集或袋类型，这个类型允许集合中的元素出现多次，是不是很神奇，一起来看一下。
from collections import Counter
inventory = Counter()
fruits = {'apple': 3, 'orange': 2}
inventory.update(fruits)
fruits_add = {'peach': 5, 'apple': 2}
inventory.update(fruits_add)
print(inventory) # Counter({'apple': 5, 'peach': 5, 'orange': 2}) 进行了累




# 类：当你的数据结构更复杂的时候，你就会构造类来封装你的数据结构，Python中用类来封装有很多方法
# 基本类：
class car:
    def __init__(self, color, mileage, automatic, length):
        self.color = color
        self.mileage = mileage
        self.automatic = automatic
        self.__length = length # 私有属性：限制对属性的直接
    @classmethod # 类方法：使用@classmethod装饰器，第一个参数是cls，表示类本身。
    def to_meter(self):
        out = self.mileage * 1000
        return out
    def producer(self, country='China'):
        return country
    def car_length(self):
        return self.__length
    @staticmethod # 静态方法：使用@staticmethod装饰器，不自动接收类或实例作为第一个参数。
    def price(RMB):
        return RMB
x = car('red', 123, True, 100) # 实例化
x.color, x.mileage, x.to_meter(), x.producer(), x.producer('UK'), x.car_length(), x.price(100) # 可访问
x.length, x.price() # 不可访问

# dataclass：数据类 是在Python 3.7才开始有的。是一个新的模块专门用于数据存储的类，使用起来非常方便，可以让你的对象实例将立即获得一些有用的功能，从而节省了一些键入和手动实现的工作：
# 定义实例变量的语法较短，无需实现该.__init__()方法。
# 数据类的实例通过自动生成的.__repr__()方法自动获得漂亮的字符串表示形式。
# 实例变量接受类型注释，从而使数据类在一定程度上能够自我记录。
from dataclasses import dataclass
@dataclass
class car:
    color:str
    mileage:float
    automatic:bool
x = car('red', 123, True)
x.color, x.mileage, x.automatic

# namedtuple 命名元组
# nametuple是Python 里面短小精悍的类，它是collection模块里面的一个库函数，对于封装短小的类非常有用。让我们可以像访问属性一样去访问类的内容。使用的效率更高更接近底层的数组的效率，而且代码易读。
from collections import namedtuple
car = namedtuple('car', 'color mileage automatic') # 定义类
x = car('red', 123, True) # 实例化
x.color, x.mileage, x.automatic # 可访问

# struct 序列化的C结构。c/c++中就有struct，在那里struct叫做结构体。在Python中也使用struct，Python正是使用struct模块执行Python值和C结构体之间的转换，从而形成Python字节对象。
# x：垫字节（不存储数据） ，b：有符号字符，B：无符号字符，
# h：有符号短整型（16位），H：无符号短整型（16位），i：有符号整型（32位），I：无符号整型（32位），
# l：有符号长整型（32位），L：无符号长整型（32位），q：有符号长整型（64位），Q：无符号长整型（64位），
# f：浮点型，d：双精度浮点型，s：字符串（使用null终止符），p：字符串（使用null终止符，与s相同），
# P：指针（一个长整型），?：布尔值
from struct import Struct
struct = Struct('i2?f') # 创建结构体
x = struct.pack(23, False, True, 100.9) # 类似于实例化，实际上把这些写成二进制的数据了
struct.unpack(x) # 恢复数据



# numpy
# NumPy 是 Python 中用于科学计算的核心库，它提供了一个高性能的多维数组对象 ndarray 和用于操作这些数组的工具和函数。

# ndarray 是 NumPy 中最基本也是最核心的数据结构，它是一个 n 维数组对象，提供了高效的大型数据集支持。
# ndarray 中的所有元素必须是相同类型的。一旦创建，ndarray 的形状（维度和大小）就不能改变。ndarray 中的元素通过基于 0 的索引访问。支持通过切片语法访问子数组。
import numpy as np
array = np.array([1, 2, 3, 4, 5]) # 通过列表创建 ndarray
matrix = np.array([[1, 2, 3], [4, 5, 6]]) # 创建多维数组

# 切片和索引
print(array[0])  # 索引，输出第一个元素，1
print(array[1:3])  # 切片，输出第二个到第三个元素，[2 3]
print(array[::2])  # 步进切片，输出每隔一个元素，[1 3 5]
print(array[array > 2])  # 布尔索引，输出大于 2 的元素，[3 4 5]

# 视图和副本 NumPy 数组的切片操作返回的是原数组的视图（view），而不是副本。这意味着，对视图的修改会影响原数组。
# 创建一个视图
view = array[1:3] # array([2, 3])
view[0] = 10
# 修改视图也修改了原数组
print(array)  # 输出：[ 1 10  3  4  5]

# 形状操作 NumPy 提供了多种函数来改变数组的形状，而不改变其数据。
reshaped = array.reshape(5, 1)  # 改变形状，将一维数组变为 5x1 的二维数组
flattened = array.flatten()  # 展平数组，将多维数组变为一维数组

# 广播（Broadcasting）NumPy 的广播机制允许不同形状的数组在数学运算中协同工作，基于一定的规则自动扩展数组的形状。
# 数组和标量的运算
result = array + 5  # 将 5 广播到 array 的形状
# 一维数组和二维数组的运算
row = np.array([1, 2, 3])
matrix += row  # 将 row 广播到 matrix 的每一行

# 唯一内存存储 NumPy 数组在内存中是连续存储的，这使得它们可以利用现代 CPU 的向量化指令进行快速计算。
# 内存映射文件 NumPy 还支持内存映射文件（memory-mapped files），这允许数组直接映射到磁盘文件的一段区域，从而处理大于内存大小的大型数据集。
# 创建一个内存映射文件
mmap = np.memmap('mmap.dat', dtype='float32', mode='w+', shape=(1000, 1000))




# pandas 
# Pandas 是 Python 中用于数据分析和操作的一个强大库，它提供了多种数据结构，主要用于处理表格数据。Pandas 的主要数据结构包括 Series 和 DataFrame。
# Series 是一个一维数组结构，可以存储任何类型的数据（整数、字符串、浮点数、Python 对象等），且每个元素都有一个标签（即索引）。
# 索引：Series 中的每个元素都有一个索引，索引标签可以是数字、字符串或日期。
# 数据类型：Series 中的所有元素都是相同类型的，但 Series 本身可以包含不同类型的数据。
# 大小不变性：Series 的大小（长度）是固定的，不能改变。

# 创建 Series
import pandas as pd
s = pd.Series([1, 2, 3, 4, 5]) # 通过列表创建 Series
s_dict = pd.Series({'a': 1, 'b': 2, 'c': 3}) # 通过字典创建 Series，字典的键成为索引

# 访问 Series 数据
# 通过索引访问
print(s[0])  # 输出：1
# 通过索引标签访问
print(s_dict['a'])  # 输出：1

# 数学运算 Series 支持大量的数学运算方法，可以对整个序列进行计算或者与标量进行运算：
s = pd.Series([1, 2, 3, 4, 5])
result = s + 5  # 对每个元素加5
result = s * 2  # 对每个元素乘以2
result = s ** 2  # 对每个元素求平方

# 统计运算 Series 提供了多种统计方法，用于计算描述性统计数据：
sum_val = s.sum()
mean_val = s.mean()
median_val = s.median()
max_val = s.max()
min_val = s.min()
std_dev = s.std()
var = s.var()

# 布尔操作 Series 支持布尔索引，可以基于条件过滤数据：
# 条件过滤
filtered_s = s[s > 2]  # 返回大于2的元素

# 字符串操作 当 Series 中的元素是字符串时，Pandas 提供了一系列字符串操作方法：
s_str = pd.Series(['hello', 'world', 'pandas'])
concatenated = s_str.str.cat(sep=' ')  # 字符串连接
upper_case = s_str.str.upper() # 大小写转换
contains_p = s_str.str.contains('p') # 字符串包含检查

# 缺失数据处理 Pandas 提供了处理缺失数据的方法：
s_missing = pd.Series([1, 2, None, 4])
filled_s = s_missing.fillna(0) # 填充缺失值
dropped_s = s_missing.dropna() # 删除缺失值

# 重索引和对齐 Series 可以被重索引，以改变其索引标签：
s_reindexed = s.reindex([0, 1, 2, 3, 4, 5])

# 累积和移动窗口 Series 提供了累积和移动窗口功能，用于计算累积和或滑动窗口统计量：
cumsum = s.cumsum() # 累积和
moving_average = s.rolling(window=2).mean() # 移动平均

# 唯一值和成员资格测试 Series 提供了检查唯一值和测试成员资格的方法：
unique_values = s.unique() # 唯一值
is_in_s = 3 in s # 值是否在 Series 中

# 值计数 Series 可以计算每个值出现的次数：
value_counts = s.value_counts()

# 时间序列功能 如果 Series 的索引是日期时间类型，那么它还支持许多时间序列功能：
s_time = pd.Series([1, 2, 3], index=pd.date_range('20210101', periods=3))
s_time = s_time.asfreq('D') # 设置时间频率
s_time_shifted = s_time.shift(1) # 时间位移

# DataFrame 
# DataFrame 是一个二维表格型数据结构，可以被看作是由多个 Series 组成的（每个 Series 作为 DataFrame 的一列），且所有 Series 共享索引。
# 索引和列：DataFrame 有行索引和列名，可以快速访问数据。
# 异构类型：DataFrame 的每一列可以是不同的数据类型。
# 灵活的形状：DataFrame 的行和列可以很容易地添加或删除。

# 创建 DataFrame
# 通过字典创建 DataFrame，字典的键成为列名
df = pd.DataFrame({
    'Column1': [1, 2, 3],
    'Column2': ['a', 'b', 'c']
})
# 通过列表的列表创建 DataFrame
df_list = pd.DataFrame([
    [1, 'a'],
    [2, 'b'],
    [3, 'c']
])

# 访问 DataFrame 数据
# 通过列名访问列
print(df['Column1'])  # Series([1, 2, 3])
# 使用 .loc[] 和 .iloc[] 访问行和列
print(df.loc[0])  # 访问第一行，返回一个 Series
print(df.iloc[0])  # 也访问第一行，返回一个 Series

# 操作 DataFrame Pandas 提供了丰富的方法来操作 DataFrame，包括：
# 选择和过滤：根据条件选择和过滤数据。
# 聚合：对数据进行聚合操作，如求和、平均、最大/最小值等。
# 合并和连接：将多个 DataFrame 或 Series 合并成一个。
# 分组：对数据进行分组，并对每组应用函数。
# 重塑：改变 DataFrame 的形状，如透视表和交叉表。
filtered = df[df['Column1'] > 1] # 选择和过滤
sum = df['Column1'].sum() # 聚合
df_concat = pd.concat([df, df_list], axis=1) # 合并
grouped = df.groupby('Column1') # 分组
pivot_table = df.pivot_table(index='Column1', columns='Column2', values='Column1', aggfunc='sum') # 重塑
# 索引对象 Pandas 还提供了 Index 对象，它是 Series 和 DataFrame 索引的基础。
# Index 对象支持高效的数据检索，并提供了丰富的方法来操作索引。

# 数据清洗和预处理
df.dropna() # 删除缺失值
df.fillna(value) # 填充缺失值
df.replace(to_replace, value) # 替换值
df.astype(type) # 数据类型转换：

# 数据选择和过滤
df[df['Column'] > value] # 条件筛选
df.query('Column > @value') # 使用 query 方法
df.eval('Column1 + Column2 > @value') # 使用 eval 表达式

# 数据聚合和摘要
df.describe() # 行或列的描述性统计
df.groupby('Column').agg(['sum', 'mean', 'max']) # 聚合函数
pd.pivot_table(df, values='Column', index='Column1', columns='Column2', aggfunc='sum') # 透视表

# 数据合并和连接
pd.concat([df1, df2], axis=1) # concat 合并
pd.merge(df1, df2, on='Key', how='inner') # merge 连接
df1.join(df2, on='Key') # join 连接

# 数据重塑
df.melt(id_vars, value_vars, var_name, value_name) # melt 函数
df.pivot(index, columns, values) # pivot 函数

# 字符串操作 
df['Column'].str.contains(pattern) # 对 DataFrame 中的字符串列进行操作
df['Column'].str.cat(sep=' ')

# 时间序列功能
# 如果 DataFrame 包含日期时间类型的列，可以进行时间序列分析：
df.set_index('DateTimeColumn')
df.resample('D').mean()
df.shift(periods)

# 排序和排名
df.sort_values(by='Column', ascending=True) # 排序
df.rank(method='average') # 排名

# 唯一值和成员资格测试
df['Column'].unique() # 唯一值
'Value' in df['Column'] # 成员资格测试

# 值计数
df['Column'].value_counts()

# 应用函数 对 DataFrame 的列或行应用自定义函数：
df.apply(func)
df.apply(func, axis=1)  # 对行应用

# 处理大数据集
pd.read_csv('large_dataset.csv', chunksize=1000) # 分块处理

# 优化内存使用：
df = df.astype({'Column': 'category'})