#### 字符串 ####
# https://www.runoob.com/python3/python3-string.html

# 访问字符串中的值
# Python 不支持单字符类型，单字符在 Python 中也是作为一个字符串使用。
# Python 访问子字符串，可以使用方括号 [] 来截取字符串，字符串的截取的语法格式如下：
# 变量[头下标:尾下标]
# 索引值以 0 为开始值，-1 为从末尾的开始位置。

# Python 转义字符
# \(在行尾时) 续行符	
print("line1 \
    ... line2")
print("\\") # \\ 反斜杠符号	
print('\'') # \' 单引号	
print("\"") # \" 双引号	
print("\a") # \a 响铃，执行后电脑有响声
print("Hello \b World!") # \b 退格(Backspace)	
print("\000") # \000 空	
print("\n") # \n 换行
print("Hello \v World!") # \v 纵向制表符
print("Hello \t World!") # \t 横向制表符	
print('google runoob taobao\r123456') # \r 回车，将 \r 后面的内容移到字符串开头，并逐一替换开头部分的字符，直至将 \r 后面的内容完全替换完成。
print("Hello \f World!") # \f 换页	
print("\110\145\154\154\157\40\127\157\162\154\144\41") # \yyy 八进制数，y 代表 0~7 的字符，例如：\012 代表换行。	
print("\x48\x65\x6c\x6c\x6f\x20\x57\x6f\x72\x6c\x64\x21") # \xyy 十六进制数，以 \x 开头，y 代表的字符，例如：\x0a 代表换行	
# \other 其它的字符以普通格式输出

# 使用 \r 实现百分比进度：\r的作用是将光标移动到当前行的开头，而不换行
import time
for i in range(101):
    print("\r{:3}%".format(i), end='') #{:3}指定宽度为3的字符
    time.sleep(0.05)


# 字符串运算符
# + # 字符串连接
# * # 重复输出字符串
# [] # 通过索引获取字符串中字符
# [ : ] # 截取字符串中的一部分，遵循左闭右开原则，str[0:2] 是不包含第 3 个字符的
# in	# 成员运算符 - 如果字符串中包含给定的字符返回 True
# not in # 成员运算符 - 如果字符串中不包含给定的字符返回 True
# r/R # 原始字符串 - 原始字符串：所有的字符串都是直接按照字面的意思来使用，没有转义特殊或不能打印的字符。 原始字符串除在字符串的第一个引号前加上字母 r（可以大小写）以外，与普通字符串有着几乎完全相同的语法。	
# % # 格式字符串	请看下一节内容。


# Python 字符串格式化
# Python 支持格式化字符串的输出 。尽管这样可能会用到非常复杂的表达式，但最基本的用法是将一个值插入到一个有字符串格式符 %s 的字符串中。
# 在 Python 中，字符串格式化使用与 C 中 sprintf 函数一样的语法。
print ("我叫 %s 今年 %d 岁!" % ('小明', 10))

# python字符串格式化符号:
# %c # 格式化字符及其ASCII码
# %s # 格式化字符串
# %d # 格式化整数
# %u # 格式化无符号整型
# %o # 格式化无符号八进制数
# %x # 格式化无符号十六进制数
# %X # 格式化无符号十六进制数（大写）
# %f # 格式化浮点数字，可指定小数点后的精度
# %e # 用科学计数法格式化浮点数
# %E # 作用同%e，用科学计数法格式化浮点数
# %g # %f和%e的简写
# %G # %f 和 %E 的简写
# %p # 用十六进制数格式化变量的地

# 格式化操作符辅助指令:
# * # 定义宽度或者小数点精度
# - # 用做左对齐
# + # 在正数前面显示加号( + )
# <sp> # 在正数前面显示空格
# # # 在八进制数前面显示零('0')，在十六进制前面显示'0x'或者'0X'(取决于用的是'x'还是'X')
# 0 # 显示的数字前面填充'0'而不是默认的空格
# % # '%%'输出一个单一的'%'
# (var) # 映射变量(字典参数)
# m.n. # m 是显示的最小总宽度,n 是小数点后的位数(如果可用的话)
# Python2.6 开始，新增了一种格式化字符串的函数 str.format()，它增强了字符串格式化的功能。

# f-string
# f-string 是 python3.6 之后版本添加的，称之为字面量格式化字符串，是新的格式化字符串的语法。之前我们习惯用百分号 (%):
x = 'name' ; print('%s' % x)
# f-string 格式化字符串以 f 开头，后面跟着字符串，字符串中的表达式用大括号 {} 包起来，它会将变量或表达式计算后的值替换进去
x = 'name' ; print(f'{x}') 
print(f'{1 + 2}') # 表达式
print(f'{x + "2"}') # 表达式
# 在 Python 3.8 的版本中可以使用 = 符号来拼接运算表达式与结果：
print(f'{x + "2"=}')


# 字符串内建函数
capitalize() # 将字符串的第一个字符转换为大写
center(width, fillchar) # 返回一个指定的宽度 width 居中的字符串，fillchar 为填充的字符，默认为空格。
count(str, beg= 0,end=len(string)) # 返回 str 在 string 里面出现的次数，如果 beg 或者 end 指定则返回指定范围内 str 出现的次数
bytes.decode(encoding="utf-8", errors="strict") # Python3 中没有 decode 方法，但我们可以使用 bytes 对象的 decode() 方法来解码给定的 bytes 对象，这个 bytes 对象可以由 str.encode() 来编码返回。
encode(encoding='UTF-8',errors='strict') # 以 encoding 指定的编码格式编码字符串，如果出错默认报一个ValueError 的异常，除非 errors 指定的是'ignore'或者'replace'
endswith(suffix, beg=0, end=len(string)) # 检查字符串是否以 suffix 结束，如果 beg 或者 end 指定则检查指定的范围内是否以 suffix 结束，如果是，返回 True,否则返回 False。
expandtabs(tabsize=8) # 把字符串 string 中的 tab 符号转为空格，tab 符号默认的空格数是 8 。
find(str, beg=0, end=len(string)) # 检测 str 是否包含在字符串中，如果指定范围 beg 和 end ，则检查是否包含在指定范围内，如果包含返回开始的索引值，否则返回-1
index(str, beg=0, end=len(string)) # 跟find()方法一样，只不过如果str不在字符串中会报一个异常。
isalnum() # 检查字符串是否由字母和数字组成，即字符串中的所有字符都是字母或数字。如果字符串至少有一个字符，并且所有字符都是字母或数字，则返回 True；否则返回 False。
isalpha() # 如果字符串至少有一个字符并且所有字符都是字母或中文字则返回 True, 否则返回 False
isdigit() # 如果字符串只包含数字则返回 True 否则返回 False
islower() # 如果字符串中包含至少一个区分大小写的字符，并且所有这些(区分大小写的)字符都是小写，则返回 True，否则返回 False
isnumeric() # 如果字符串中只包含数字字符，则返回 True，否则返回 False
isspace() # 如果字符串中只包含空白，则返回 True，否则返回 False
istitle() # 如果字符串是标题化的(见 title())则返回 True，否则返回 False
isupper() # 如果字符串中包含至少一个区分大小写的字符，并且所有这些(区分大小写的)字符都是大写，则返回 True，否则返回 False
join(seq) # 以指定字符串作为分隔符，将 seq 中所有的元素(的字符串表示)合并为一个新的字符串
len(string) # 返回字符串长度
# ljust(width [, fillchar]) # 返回一个原字符串左对齐,并使用 fillchar 填充至长度 width 的新字符串，fillchar 默认为空格。 x.ljust(10, '#')
lower() # 转换字符串中所有大写字符为小写.
lstrip() # 截掉字符串左边的空格或指定字符。
maketrans() # 创建字符映射的转换表，对于接受两个参数的最简单的调用方式，第一个参数是字符串，表示需要转换的字符，第二个参数也是字符串表示转换的目标。
translate(table, deletechars="") # 根据 table 给出的表(包含 256 个字符)转换 string 的字符, 要过滤掉的字符放到 deletechars 参数中
x = 'name' ; trans = str.maketrans('na', 'NA') ; x.translate(trans) # NAme
max(str) # 返回字符串 str 中最大的字母。
min(str) # 返回字符串 str 中最小的字母。
# replace(old, new [, max]) # 把 将字符串中的 old 替换成 new,如果 max 指定，则替换不超过 max 次。[, max] 表示第三个参数是max，但是是可选的，不一定必须得选
rfind(str, beg=0,end=len(string)) # 类似于 find()函数，不过是从右边开始查找
rindex( str, beg=0, end=len(string)) # 类似于 index()，不过是从右边开始
# rjust(width, [, fillchar]) # 返回一个原字符串右对齐,并使用fillchar(默认空格）填充至长度 width 的新字符串
rstrip() # 删除字符串末尾的空格或指定字符
split(str="", num=string.count(str)) # 以 str 为分隔符截取字符串，如果 num 有指定值，则仅截取 num+1 个子字符串
splitlines([keepends]) # 按照行('\r', '\r\n', \n')分隔，返回一个包含各行作为元素的列表，如果参数 keepends 为 False，不包含换行符，如果为 True，则保留换行符。
startswith(substr, beg=0,end=len(string)) # 检查字符串是否是以指定子字符串 substr 开头，是则返回 True，否则返回 False。如果beg 和 end 指定值，则在指定范围内检查。
strip([chars]) # 在字符串上执行 lstrip()和 rstrip()
swapcase() # 将字符串中大写转换为小写，小写转换为大写
title() # 返回"标题化"的字符串,就是说所有单词都是以大写开始，其余字母均为小写(见 istitle())
upper() # 转换字符串中的小写字母为大写
zfill(width) # 返回长度为 width 的字符串，原字符串右对齐，前面填充0
isdecimal() # 检查字符串是否只包含十进制字符，如果是返回 true，否则返回 false