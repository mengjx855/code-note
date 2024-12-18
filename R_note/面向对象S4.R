# Jinxmeng, 20240809, 20240809 -------------

setwd("/R_proj/Jinxin/S4class/")

# 一些程序员认为S3类不具有面向对象编程固有的安全性。例如，你可以任意修改S3类，哪怕是不合法的修改。相比而言，S4类更加安全。

# S4类的创建
# 可以调用setClass()来定义一个S4类。
setClass(Class, representation, prototype, contains=character(),
         validity, access, where, version, sealed, package,
         S3methods = FALSE, slots)

# Class: 定义类名
# slots: 定义属性和属性类型
# prototype: 定义属性的默认值
# contains=character(): 定义父类，继承关系
# validity: 定义属性的类型检查
# where: 定义存储空间
# sealed: 如果设置TRUE，则同名类不能被再次定义
# package: 定义所属的包

# 创建一个S4对象实例
# 为了方便我们检查对象的类型，引入pryr包作为辅助工具。
library(pryr)

# 定义一个S4对象
setClass("Person",slots=list(name="character",age="numeric"))
# 实例化一个Person对象
father<-new("Person",name="F",age=44)
# 查看father对象，有两个属性name和age
father
# 查看father对象类型，为Person
class(father)

# 创建一个有继承关系的S4对象
# 创建一个S4对象Person
setClass("Person",slots=list(name="character",age="numeric"))
# 创建Person的子类
setClass("Son",slots=list(father="Person",mother="Person"),contains="Person")
# 实例化Person对象
father<-new("Person",name="F",age=44)
mother<-new("Person",name="M",age=39)
# 实例化一个Son对象
son<-new("Son",name="S",age=16,father=father,mother=mother)
# 查看son对象的name属性
  son@name
# 查看son对象的age属性
son@age
# 查看son对象的father属性
son@father
# 查看son对象的mother属性
slot(son,"mother")
# 检查son类型
otype(son)
# 检查son@name属性类型
otype(son@name)
# 检查son@mother属性类型
otype(son@mother)
# 用isS4()，检查S4对象的类型
isS4(son)
isS4(son@name)
isS4(son@mother)

# S4对象的默认值
setClass("Person",slots=list(name="character",age="numeric"))
# 属性age为空
a<-new("Person",name="a")
a
# 设置属性age的默认值20
setClass("Person",slots=list(name="character",age="numeric"),prototype = list(age = 20))
# 属性age为空
b<-new("Person",name="b")
# 属性age的默认值是20
b

# S4对象的类型检查
setClass("Person",slots=list(name="character",age="numeric"))
# 传入错误的age类型
bad<-new("Person",name="bad",age=3)
# 设置age的非负检查
setValidity("Person",function(object) {
  if (object@age <= 0) stop("Age is negative.")
  }
)
# 修传入小于0的年龄
bad2<-new("Person",name="bad",age=-1)

# 从一个已经实例化的对象中创建新对象
# S4对象，还支持从一个已经实例化的对象中创建新对象，创建时可以覆盖旧对象的值
setClass("Person",slots=list(name="character",age="numeric"))
# 创建一个对象实例n1
n1<-new("Person",name="n1",age=19);n1
# 从实例n1中，创建实例n2，并修改name的属性值
n2<-initialize(n1,name="n2");n2

# 访问对象的属性
# 在S3对象中，一般我使用$来访问一个对象的属性，但在S4对象中，我们只能使用@来访问一个对象的属性
setClass("Person",slots=list(name="character",age="numeric"))
a<-new("Person",name="a")
# 访问S4对象的属性
a@name
slot(a, "name")
# 错误的属性访问
a$name
a[1]


# S4的泛型函数
# 在R语言中，S4类是一种更为严格和正式的面向对象编程机制，与之相对的是较为简单和灵活的S3类。S4类允许用户定义更复杂的数据结构和方法，提供了更强的类型检查和验证机制。
# 泛型函数（Generic Function）
# 泛型函数在S4系统中有着重要的地位。
# 泛型函数本质上是一个占位符，它定义了一个函数接口，但不具体实现该函数的行为。
# 其行为由一个或多个具体方法（Methods）来实现。
# 具体的方法是根据输入数据的类（class）来选择的。

# S4的泛型函数实现有别于S3的实现，S4分离了方法的定义和实现，如在其他语言中我们常说的接口和实现分离。
# 通过setGeneric()来定义接口，通过setMethod()来定义现实类。这样可以让S4对象系统，更符合面向对象的特征。

# 普通函数的定义和调用
work<-function(x) cat(x, "is working")
work('Conan')

# 定义Person对象
setClass("Person",slots=list(name="character",age="numeric"))
# 定义泛型函数work，即接口
setGeneric("work",function(object) standardGeneric("work"))

# 定义work的现实，并指定参数类型为Person对象
setMethod("work", signature(object = "Person"), function(object) cat(object@name , "is working") )
# 创建一个Person对象a
a<-new("Person",name="Conan",age=16)
# 把对象a传入work函数
work(a)

# 通过S4对象系统，把原来的函数定义和调用2步，为成了4步进行：
# （1）定义数据对象类型 class
# （2）定义接口函数 setGeneric
# （3）定义实现函数 setMethod
# （4）把数据对象以参数传入到接口函数，执行实现函数
# 通过S4对象系统，是一个结构化的，完整的面向对象实现。

# 查看S4对象的函数
# 当我们使用S4对象进行面向对象封装后，我们还需要能查看到S4对象的定义和函数定义。
# 还以上面Person和work的例子
# 检查work的类型
ftype(work)
# 直接查看work函数
work
# 查看work函数的现实定义
showMethods(work)
# 查看Person对象的work函数现实
getMethod("work", "Person")
selectMethod("work", "Person")

# 检查Person对象有没有work函数
existsMethod("work", "Person")
hasMethod("work", "Person")
