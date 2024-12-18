#### os ####
# https://zhuanlan.zhihu.com/p/150835193
# https://www.runoob.com/python/os-file-methods.html


# os就是“operating system”的缩写，顾名思义，os模块提供的就是各种 Python 程序与操作系统进行交互的接口。
# 通过使用os模块，一方面可以方便地与操作系统进行交互，另一方面页可以极大增强代码的可移植性。如果该模块中相关功能出错，会抛出OSError异常或其子类异常。
# 注意，如果是读写文件的话，建议使用内置函数open()；
# 如果是路径相关的操作，建议使用os的子模块os.path；
# 如果要逐行读取多个文件，建议使用fileinput模块；
# 要创建临时文件或路径，建议使用tempfile模块；
# 要进行更高级的文件和路径操作则应当使用shutil模块。
# 当然，使用os模块可以写出操作系统无关的代码并不意味着os无法调用一些特定系统的扩展功能，但要切记一点：一旦这样做就会极大损害代码的可移植性。
# 此外，导入os模块时还要小心一点，千万不要为了图调用省事儿而将os模块解包导入，也就是不要使用from os import *来导入os模块；否则os.open()将会覆盖内置函数open()，从而造成预料之外的错误。

# 注意，os模块中大多数接受路径作为参数的函数也可以接受“文件描述符”作为参数。
# 文件描述符：file descriptor，在 Python 文档中简记为 fd，文件描述符是一个整数，它唯一地标识了系统中打开的文件。当你打开一个文件时，操作系统会返回一个文件描述符，这个描述符用于后续对该文件的所有操作，直到文件被关闭。
# 也就类似于perl中的文件句柄。
fd = os.open(path, os.O_RDONLY) # 打开文件并获取文件描述符


# os.name
# 该属性宽泛地指明了当前 Python 运行所在的环境，实际上是导入的操作系统相关模块的名称。这个名称也决定了模块中哪些功能是可用的，哪些是没有相应实现的。
# 目前有效名称为以下三个：posix，nt，java。
# 其中posix是 Portable Operating System Interface of UNIX（可移植操作系统接口）的缩写。Linux 和 Mac OS 均会返回该值；
# nt全称应为“Microsoft Windows NT”，大体可以等同于 Windows 操作系统，因此 Windows 环境下会返回该值；
# java则是 Java 虚拟机环境下的返回值。
os.name
sys.platform # 可以返回更详细的信息


# os.environ
# 返回环境变量的信息，是一个字典，键是环境名称，值是路径
os.environ["HOME"]


# os.walk()
# 这个函数需要传入一个路径作为top参数，函数的作用是在以top为根节点的目录树中游走，对树中的每个目录生成一个由(dirpath, dirnames, filenames)三项组成的三元组。
# 其中，dirpath是一个指示这个目录路径的字符串，dirnames是一个dirpath下子目录名（除去“.”和“..”）组成的列表，filenames则是由dirpath下所有非目录的文件名组成的列表。要注意的是，这些名称并不包含所在路径本身，要获取dirpath下某个文件或路径从top目录开始的完整路径，需要使用os.path.join(dirpath, name)。
os.system('tree xx')
# xx
# ├── xx.txt
# └── xx2
#     └── xx.txt
# 1 directory, 2 files
for i in os.walk('xx'):
    print(i)
# ('xx', ['xx2'], ['xx.txt'])
# ('xx/xx2', [], ['xx.txt'])
# 当前路径，当前路径下所有的路径（不包括'.'和'..'），当前路径下所有的文件；


os.listdir() # 列出（当前）目录下的全部路径（及文件）。该函数存在一个参数，用以指定要列出子目录的路径，默认为“.”，即“当前路径”。函数返回值是一个列表，其中各元素均为字符串，分别是各路径名和文件名。
os.mkdir() # 用处是“新建一个路径”。需要传入一个类路径参数用以指定新建路径的位置和名称，如果指定路径已存在，则会抛出FileExistsError异常。该函数只能在已有的路径下新建一级路径，否则（即新建多级路径）会抛出FileNotFoundError异常。
os.makedirs() # 在需要新建多级路径的场景下，可以使用该函数。函数os.makedirs()执行的是递归创建，若有必要，会分别新建指定路径经过的中间路径，直到最后创建出末端的“叶子路径”。
os.remove() # 用于删除文件，如果指定路径是目录而非文件的话，就会抛出IsADirectoryError异常。
os.rmdir() # 删除目录应该使用这个函数。
os.removedirs() # 递归删除目录的函数os.removedirs()，该函数会尝试从最下级目录开始，逐级删除指定的路径，几乎就是一个os.makedirs()的逆过程；一旦遇到非空目录即停止。
os.rename() # 该函数的作用是将文件或路径重命名，一般调用格式为os.rename(src, dst)，即将src指向的文件或路径重命名为dst指定的名称。
os.renames() # 和上两个函数一样，该函数也有对应的递归版本os.renames()，能够创建缺失的中间路径。注意，这两种情况下，如果函数执行成功，都会调用os.removedir()函数来递归删除源路径的最下级目录。
os.getcwd() # “getcwd”实际上是“get the current working directory”的简写，顾名思义，也就是说这个函数的作用是“获取当前工作路径”。在程序运行的过程中，无论物理上程序在实际存储空间的什么地方，“当前工作路径”即可认为是程序所在路径；与之相关的“相对路径”、“同目录下模块导入”等相关的操作均以“当前工作路径”为准。
os.chdir() # “chdir”其实是“change the directory”的简写，因此os.chdir()的用处实际上是切换当前工作路径为指定路径。其中“指定路径”需要作为参数传入函数os.chdir()，该参数既可以是文本或字节型字符串，也可以是一个文件描述符，还可以是一个广义的类路径（path-like）对象。若指定路径不存在，则会抛出FileNotFoundError异常。
os.chmod(path, mode) # 更改权限
os.chown(path, uid, gid) # 更改文件所有者
os.close(fd) # 关闭文件描述符 fd
os.open(file, flags) # 打开一个文件，并且设置需要的打开选项，mode参数是可选的
os.read(fd, n) # 从文件描述符 fd 中读取最多 n 个字节，返回包含读取字节的字符串，文件描述符 fd对应文件已达到结尾, 返回一个空字符串。