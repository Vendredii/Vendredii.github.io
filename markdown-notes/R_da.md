# 判别分析
Discriminant analysis是一种分类的手段，就是在已经确定分类的情况下根据新的值推断新的东西属于哪一类。
[TOC]
## 线性判别分析lda
使用的数据来自[Guerin Chloe](http://doi.org/10.1093/aob/mcz175)关于利用棕榈科内植物不饱和脂肪酸的各种含量推断其所属的亚族的案例，可以在链接中的supplementry中的表S7中获得。
准备数据与环境
```r
library(MASS)
seed <- read.csv("seed.csv", head = T)
head(seed)
#        class C8.0 C10.0 C12.0 C14.0 C16.0 C16.1 C17.0 C18.0
# 1 Attaleinae  9.4   8.2  49.2  13.9   6.0     0     0   3.2
# 2 Attaleinae  8.9   5.9  41.0  14.6   8.9     0     0   5.3
# 3 Attaleinae  6.0   5.2  53.8  17.9   6.3     0     0   1.7
# 4 Attaleinae 13.5  15.5  36.2   7.2   4.8     0     0   2.1
# 5 Attaleinae 15.0  13.0  44.5   8.0   3.9     0     0   1.7
# 6 Attaleinae  9.2   6.7  49.1  18.4   8.5     0     0   2.8
#   C18.1n7 C18.1n9 C18.2 C18.3n3 C18.3n6 C20.0 C20.1n9 C22.0
# 1       0     8.6   1.5       0       0   0.0       0     0
# 2       0    12.8   2.4       0       0   0.1       0     0
# 3       0     7.3   1.8       0       0   0.0       0     0
# 4       0    16.6   4.1       0       0   0.0       0     0
# 5       0     9.6   4.2       0       0   0.1       0     0
# 6       0     4.4   0.9       0       0   0.0       0     0
#   C22.1 C24.0
# 1     0     0
# 2     0     0
# 3     0     0
# 4     0     0
# 5     0     0
# 6     0     0
```
进行lda分析
```r
ord <- lda(class ~., seed)
ord
# Call:
# lda(class ~ ., data = seed)

# 各个分类数据占总体的比重，用ord$prior调用
# Prior probabilities of groups:
#        Attaleinae       Bactridinae        Dypsidinae 
#         0.2682927         0.2926829         0.1219512 
#      Livistoninae Ptychospermatinae 
#         0.1951220         0.1219512 

# 各个分类的均值向量，用ord$means调用
# Group means:
#                        C8.0    C10.0    C12.0    C14.0
# Attaleinae        10.354545 9.645455 44.27273 12.67273
# Bactridinae        2.358333 2.091667 52.99167 23.90833
# Dypsidinae         1.060000 1.360000 41.08000 21.52000
# Livistoninae       0.400000 0.512500 25.46250 12.51250
# Ptychospermatinae  0.100000 0.120000  7.40000 15.58000
#                       C16.0      C16.1  C17.0    C18.0
# Attaleinae         6.263636 0.00000000 0.0000 2.727273
# Bactridinae        6.083333 0.09166667 0.0000 1.941667
# Dypsidinae        10.900000 0.04000000 0.1000 2.020000
# Livistoninae      14.275000 0.17500000 0.0875 3.075000
# Ptychospermatinae 22.280000 0.18000000 0.1200 4.260000
#                   C18.1n7  C18.1n9     C18.2 C18.3n3
# Attaleinae           0.00 11.43636  2.572727  0.0000
# Bactridinae          0.00  7.75000  2.750000  0.0000
# Dypsidinae           0.00 11.22000  9.820000  0.2400
# Livistoninae         0.00 25.18750 17.475000  0.2625
# Ptychospermatinae    0.18 26.56000 22.100000  0.4200
#                   C18.3n6       C20.0     C20.1n9  C22.0
# Attaleinae           0.00 0.018181818 0.009090909 0.0000
# Bactridinae          0.00 0.008333333 0.008333333 0.0000
# Dypsidinae           0.04 0.180000000 0.120000000 0.1200
# Livistoninae         0.00 0.112500000 0.225000000 0.0375
# Ptychospermatinae    0.00 0.300000000 0.120000000 0.2000
#                    C22.1 C24.0
# Attaleinae        0.0000   0.0
# Bactridinae       0.0000   0.0
# Dypsidinae        0.0000   0.2
# Livistoninae      0.0625   0.1
# Ptychospermatinae 0.0000   0.2

# 降维后的矩阵，用ord$scaling调用
# Coefficients of linear discriminants:
#               LD1        LD2         LD3         LD4
# C8.0     7.455291  -2.314600   1.2084649   2.7055447
# C10.0    8.246015  -3.005224   0.7405566   2.8601837
# C12.0    7.699031  -2.804182   0.7473331   2.8390547
# C14.0    7.928075  -2.752640   0.8793914   2.9278632
# C16.0    8.274876  -2.662281   0.7989097   2.6498807
# C16.1    4.627508  -3.913803  -0.3161783   3.2083925
# C17.0   -9.159284 -38.262766  16.2293611   0.8829658
# C18.0    7.525268  -2.621863   1.2263565   3.2447562
# C18.1n7 17.025776 -14.368888  10.4090349   3.7960774
# C18.1n9  7.831366  -2.761313   0.8586662   2.9637938
# C18.2    7.969108  -2.650793   0.7212941   3.0479520
# C18.3n3  8.678906  -1.848086   2.8712302   1.6793984
# C18.3n6  3.298594 -14.894372   8.8791414  -1.3915681
# C20.0   28.289201  -5.148618  -3.5760683  -3.9201428
# C20.1n9  6.289320 -31.485980  13.9892531 -10.0222421
# C22.0   -9.382732  27.562065 -20.5810418   4.4250677
# C22.1   12.005084   9.744644  -2.5442515   7.6616473
# C24.0   10.673824   5.938224  -4.1534212  -2.1033277

# 降维后各个分量的权重
# Proportion of trace:
#    LD1    LD2    LD3    LD4 
# 0.8202 0.1070 0.0531 0.0197 
```
查看预测结果
```r
result <- predict(ord, seed)
table(seed$class, result$class)
  #                 Attaleinae Bactridinae Dypsidinae Livistoninae Ptychospermatinae
  # Attaleinae                 9           2          0          0          0
  # Bactridinae                1          11          0          0          0
  # Dypsidinae                 0           0          5          0          0
  # Livistoninae               0           0          0          8          0
  # Ptychospermatinae          0           0          0          0          5
```
### 使用ggord进行绘制
```r
library(ggord)
p <- ggord(ord, seed$class)
p
```
有点丑，凑合吧。
![hmap2](Rmodel/Rplot37.jpeg)