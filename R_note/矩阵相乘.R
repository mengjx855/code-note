# Jinxin Meng, 20241018, 20241018 ---------------------

KO_n <- read.delim("KOs.count.tsv", row.names = 1) %>% as.matrix()
species_ab <- readRDS("../profile/genomospecies.tpm.rds") %>% as.matrix()
KO_ab <- KO_n %*% species_ab

# KO_n
# name  g1  g2  g3  g4 ... gn
# K001 1    2   3   2
# K002 2    1   3   2
# K003 1    1   2   1
# K004 2    0   1   2
# Km

# species_ab
# name  s1  s2  s3  s4 ... sn
# g1    0.1 0.3 0.2 0.1
# g2    0.3 0.4 0.4 0.2
# g3    0.5 0.2 0.2 0.1
# g4    0.2 0.1 0.4 0.2

# 矩阵相乘 KO_n %*% species_ab
# name                 s1                            s2                                s3
# K001  sum([1,2,3,2]*[0.1,0.3,0.5,0.2]) sum([1,2,3,2]*[0.3,0.4,0.2,0.1]) sum([1,2,3,2]*[0.2*0.4*0.2*0.4])
# K002  sum([2,1,3,2]*[0.1,0.3,0.5,0.2]) 
# K003  sum([1,1,2,1]*[0.1,0.3,0.5,0.2]) 
# K004  sum([2,0,1,2]*[0.1,0.3,0.5,0.2]) 

a = matrix(sample(1:3, size = 9, replace = T), 3, 3)
a
     [,1] [,2] [,3]
[1,]    3    1    2
[2,]    2    3    2
[3,]    1    1    2

a %*% a
     [,1] [,2] [,3]
[1,]   13    8   12
[2,]   14   13   14
[3,]    7    6    8

