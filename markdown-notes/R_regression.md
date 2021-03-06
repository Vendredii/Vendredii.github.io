# 各种各样的回归分析
[TOC]
## 从lm( )开始进行简单一元线性回归
首先需要从gapminder上获得一些练习数据，gapminder是一个有着全球各种数据的公益网站，我们可以通过R包“gapminder”去下载(http://github.com/jennybc/gapminder)。
我们先加载一些必要的包和数据：
```r
library(tidyverse)
library(broom)
d <- gapminder::gapminder
```
tidyverse包是个合集，里面有dplyr实现数据整理，tidyr实现数据筛选，stringr实现字符串操作，还有ggplot2去作图。
broom包接受R中内置函数的杂乱输出（如lm和nls），并将它们转为整齐的数据帧。
第三行中我们使用双冒号在不打开包的情况下加载包中的某一功能，其语法为packagename::functionname。同时双冒号还可以在多个包下有同名函数时指定我们需要的包。
接下来我们打开d，这是一个6列1704行的表格
```r
d
```
然后我们对数据进行绘图：
```r
x = ggplot(d, aes(year,lifeExp))+
  geom_point(alpha = 0.5, position = position_jitter(width = 0.6))
x
```
可以得到这样一张图：  
![ggplot导出的简单散点图](R/Rplot1.jpeg)  

接着我们来尝试进行寿命与年份间的简单线性回归：最简单的lm模型可以表示为df <- lm(a<-b)，即为用b去拟合a，并将结果输出到df里，然后我们输入`summary(df)`就能看到拟合结果了。我们可以利用coef()函数读取回归系数，并通过tidy()函数将回归结果整理成表格。
```r
lm1 <- lm(lifeExp ~ year, data = d)
summary(lm1)
coef(lm1)
tidy(lm1)
##得到结果如下
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-39.949  -9.651   1.697  10.335  22.158 
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -585.65219   32.31396  -18.12   <2e-16 ***
#year           0.32590    0.01632   19.96   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 11.63 on 1702 degrees of freedom
#Multiple R-squared:  0.1898,	Adjusted R-squared:  0.1893 
#F-statistic: 398.6 on 1 and 1702 DF,  p-value: < 2.2e-16
```
如何理解这里面的各种玩意儿呢？我们直接看Coefficients部分：
其中Estimate为估值，Std.Error为标准误差，t value为T值，Pr为P值，一般而言，P<0.05以为数据有显著差异，我们认为这些数据通过了显著性检验，可以用了，有了统计学意义，而P<0.001意为数据有极显著差异（现在一般不提这个，误导性比较强）。
接下来可以看看Multiple R-squared（拟合优度）和Adjusted R-squared（修正的拟合优度），就是所谓的R方了，数值越高说明拟合程度越好。
而F-statistic则是F检验，是用来检验方程的整体显著性的，我们通过观察它的P值去看看方程整体是否显著。
### 关于P-value和R^2
以前一直说什么显著差异极显著差异的，这时候我tm就要问了，我们线性回归为什么要扯什么显著差异，又不是做因素分析。所以P值有什么P用呢？
P-value是拒绝原假设犯第一类假设错误的概率（？），百度说是“原假设是正确的，但我们却拒绝了原假设”（？）。举个例子就是假设抛均匀硬币正面的概率是50%（P=0.5），那么现在抛5次硬币都是正面的概率显然为0.5^5=0.03125，如果我抛5次硬币都是正面，那么P=0.03125<0.05，有统计学意义了，就可以得出结论推翻原先的关于均匀硬币的说法了，我们据此可以认为这个硬币是不均匀的。
那么在线性回归中，原假设即为“其实这些数据是随机的，根本不存在什么线性”，我们通过P值检验推翻了这一原假设，即可正面“这些数据真不是随机的，而是有关系的”。
R方是什么，在SPSS中我们可以看到R-squared = SSR/TSS = 1-RSS/TSS，其中SSR为解释方差，RSS为残差平方和，TSS为固有方差。那么显然，可以用公式来表示这个R-squared：
$$R^{2}=1-\frac{\sum_{i}\left (y^{(i)}-\hat{y}^{(i)}\right )^2}{\sum_{i}\left (y^{(i)}-\bar{y}\right )^2}$$
对于$\sum_{i}\left (y^{(i)}-\hat{y}^{(i)}\right )^2$，我们已经很熟悉了([见生态学数学原理线性代数章节](./MathPrinciples.html))，而对于$\sum_{i}\left (y^{(i)}-\bar{y}\right )^2$，这就是使用平均数来预测产生的错误（损失函数）（在ML中称为基准模型（Baseline Model），那么如果说我们辛辛苦苦回归来的错误甚至多于随便求个平均值的错误，$R^2$就会小于0，就白给了。所以R方的值应当在0-1之间，且越大说明预测越准确。
那么调整R方又是什么呢？在R方的计算中，不断增加变量会提升模型的效果，但是其实并没有什么效果，而调整R方可以惩罚那些不显著的变量，来略微调低原先的R方。

回到正题，我们继续对gapminder数据进行一些处理并回归，在tidyverse包中，我们可以通过管道符%>%将前一个命令的输出作为后一个命令的输入，而不是使用嵌套函数搞一堆简称出来：
```r {class=line-numbers}
#我们先将上面的数据d称为life（套娃），然后按国家进行分类
life <- d %>% 
  group_by(country) %>%
#接着把分类完的数据使用summarise这个统计描述函数，将lifeExp定义为算术平均寿命，将gdpPercap定义为算术平均gdp。
  summarise(lifeExp = mean(lifeExp), gdpPercap = mean(gdpPercap)) 
#对数据带了log，取了对数，然后画成散点图
ggplot(life, aes(log(gdpPercap), log(lifeExp))) + geom_point()
#回归
lm2 <- lm(log(lifeExp) ~ log(gdpPercap), data = life)
summary(lm2)
##接着系统会生成结果：
#lm(formula = log(lifeExp) ~ log(gdpPercap), data = life)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.42480 -0.05350  0.01827  0.05729  0.23078 
#
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    2.920294   0.060789   48.04   <2e-16 ***
#log(gdpPercap) 0.139064   0.007296   19.06   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1036 on 140 degrees of freedom
#Multiple R-squared:  0.7219,	Adjusted R-squared:  0.7199 
#F-statistic: 363.3 on 1 and 140 DF,  p-value: < 2.2e-16
```
![ggplot导出的简单散点图](R/Rplot2.jpeg)

我们注意一下现在的$R^2$是**0.7219**，结果很好，所以……

### 为什么要取对数？

我们先试试上面的数据不取对数会如何？
代码就不列了，结果如下：
```r
#lm(formula = lifeExp ~ gdpPercap, data = life)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-42.058  -5.709   1.995   6.082  12.468 
#
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 5.308e+01  9.106e-01   58.29   <2e-16 ***
#gdpPercap   8.862e-04  8.192e-05   10.82   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 8.255 on 140 degrees of freedom
#Multiple R-squared:  0.4553,	Adjusted R-squared:  0.4514 
#F-statistic:   117 on 1 and 140 DF,  p-value: < 2.2e-16
```
![ggplot导出的简单散点图](R/Rplot3.jpeg)

此时不仅仅是图丑出天际的问题了，而是$R^2$变成了**0.4553**，数据的拟合程度都下降了好多！
首先要说明的是对数函数在其定义域内单调增，所以取对数后不会改变数据的相对关系，因此不用担心取了对数后数据不严谨了，乱了之类的。
接下来我尝试解释一下为什么取对数后我们的R方会提高，这也许可以从百度搜索引擎的机制开始讲起：
如果要搜索我们想要的信息，比如想去搜“薛定谔的猫”，那么搜索引擎是这么找到目标的呢？我们易得如果一个词在网页中出现的次数越少那么它越重要，因为显然我们搜索出的结果中肯定是首先与“薛定谔”有关的，而不是与“的”有关的，或者与“猫”有关的。也就是说如果一个关键词只在很少的网页中出现，我们通过它就容易锁定搜索目标，它的权重也就应该大。反之如果一个词在大量网页中出现，我们看到它仍然不很清楚要找什么内容，因此它应该小。概括地讲，假定一个关键词$w$在$Dw$个网页中出现过，那么$Dw$越大，$w$的权重越小，反之亦然。
这里就引入了“逆文本频率指数”（Inverse document frequency缩写为IDF）的概念，它的公式为
$$IDF=log(D/Dw)$$
<center>（lgo就是ln）,其中D是全部网页数。</center>
假定现在有D=10亿个网页，且“的”在所有网页中都出现，那么IDF(的)=log(10亿/1o亿)=0。而“薛定谔”只在200万个网页中出现，那么IDF(薛定谔)=6.2，而我们假设猫在其中1亿个网页都出现了，所以IDF(猫)=2.3，那么我们可以认为在“薛定谔的猫”中，“薛定谔”贡献最大，“猫”也有贡献，而“的”其实没用，这也是符合我们现实逻辑的。
由此我们可以得知，取对数可以将数据在整个值域中因不同区间而带来的差异降到最小。而且可以改变变量的尺度，使得数据更加平稳。

### 函数的拟合优度：从熵到AIC

信息熵反映了一个系统的有序化程度，一个系统越是有序，那么它的信息熵就越低，反之就越高，以下为熵的定义：
如果一个随机变量$X$的可能取值为$X=\left \{x_1,x_2,...,x_n\right \}$，对应的概率为$p(X=x_i)(i=1,2,...,n)$，则随机变量X的熵为：
$$H(x)=-\sum_{i=1}^{n}p(x_i)\textup{log}p(x_i)$$
那么相对熵又称交叉熵，Kullback-Leible散度（KL散度），设$p(x)$和$q(x)$式X取值的两个概率分布，则$p$对$q$的相对熵为：
$$D(p||q)=\sum_{i=1}^{n}p(x)\textup{log}\frac{p(x)}{q(x)}$$
在一定程度上，熵可以度量两个随机变量的距离，KL散度是衡量两个概率分布$p$和$q$的非对称性度量。
该相对熵有以下两个性质：
1.KL散度并非两者的距离函数，因为它是非对称的：$D(p||q)\neq D(q||p)$
2.相对熵的值非负（可通过吉布斯不等式证明）
Akaike发现K-L距离的估计在实际情况中存在着过估计，且过估计的量近似等于需要估计的模型的参数个数K+1。于是他进行了优化，并定义了最小信息化准则AIC作为模型挑选的准则(L为模型的极大似然函数)：
$$AIC=2k-2ln(L)$$
当模型的误差服从独立正态分布时：
$$AIC=nlog(\hat{\sigma}^2)+2(k+1)$$
其中$\hat{\sigma}^2=\frac{RSS}{n}$，k为参数个数，$\hat{\sigma}^2$是$\sigma^2$的极大似然估计，$n$为样本大小，$RSS$为残差平方和。
AIC为模型选择提供了有效的规则，但也有不足之处。当样本容量很大时，在AIC准则中拟合误差提供的信息就要受到样本容量的放大，而参数个数的惩罚因子却和样本容量没关系，因此当样本容量很大时，使用AIC准则选择的模型不收敛与真实模型，它通常比真实模型所含的未知参数个数要多。BIC（Bayesian InformationCriterion）贝叶斯信息准则是Schwartz在1978年根据Bayes理论提出的判别准则，称为SBC准则(也称BIC)，弥补了AIC的不足。SBC的定义为： 
$$BIC = ln(k) - 2ln(L)$$

## 分位数回归

分位数回归可以用R包quantreg实现：
```r
library(SparseM)
library(quantreg)
data("engel")
#engel是quantreg中的自由数据，有235条，2个变量，一个是income收入，一个是foodexp食品支出，我们可以借此考察收入与食品支出的关系
#建立一个0.5分位数回归，rq函数就是quantreg里进行分位数回归的函数，tau即为分位数值
#是不是可以这么理解：所谓的α分位数回归，就是希望回归曲线之下能够包含α（一个百分数）的数据点？
fit1 <- rq(foodexp ~ income, tau = .5, data = engel)
fit1
summary(fit1)
##可得到：
#Call: rq(formula = foodexp ~ income, tau = 0.5, data = engel)
#
#tau: [1] 0.5
##coefficient列给出了估计的截距和斜率
##lower bd和upper bd则是估计的置信区间
#Coefficients:
#            coefficients lower bd  upper bd 
#(Intercept)  81.48225     53.25915 114.01156
#income        0.56018      0.48702   0.60199
```
```r {class=line-numbers}
#接着使用其内置的函数计算模型拟合残差，并绘制残差图
r1 <- resid(fit1)
plot(r1)
#得到结果表明残差均匀分布在(-200,200)区间内，说明拟合效果还不错
#绘制income和foodexp的散点图，并绘制不同分位数的分位数回归曲线：
#attach命令可以避免通过$来每次调用数据框中的变量
attach(engel)
#绘制散点图
plot(income, foodexp, cex = 0.25, type = "n", xlab ="Household Income", ylab ="Food Expenditure")
points(income,foodexp , cex = 0.25, col = "blue")
#通过abline绘制直线，这里绘制0.5分位回归线
abline(rq(foodexp ~income , tau = 0.5), col = "blue")
#绘制线性回归线
abline(lm(foodexp ~income  ), lty = 2, col = "red")
#给定下列分位数
taus <- c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
#分别分位数回归
for (i in 1:length(taus)) {
  abline(rq(foodexp ~income  , tau = taus[i]), col = "gray")
}
#接触attach()的绑定
detach(engel)
```
![分位数回归1](R/Rplot26.jpeg)

接下来探讨分位数回归与恩格尔系数的关系：
```r
#使用within函数添加一列xx变量，xx为不同人群不同收入的占比
engel<-within(engel,xx <- income - mean(income))
fit1 <- summary(rq(foodexp~xx,tau=2:98/100,data = engel))
#通过mfrow功能实现1页放2幅图
plot(fit1,mfrow = c(1:2))
```
![分位数回归2](R/Rplot27.jpeg)

上图绘制了0.02到0.98这个区间中，每隔0.01做一次分位回归，其中黑色实心点代表回归曲线的截距值，阴影部分代表95%置信区间，红色实线和虚线分别代表的是，线性回归曲线的截距值和置信区间。从图中可以看出，收入对于0.02分位的foodexp的影响在0.35左右，对于0.98分位的影响在0.7左右。即为不同的分位数回归对应着不同的置信区间。
不同分位数的分位回归的截距值是否是由于抽样误差造成的，我们同样需要假设检验进行验证，那么我们使用Wald检验进行验证。得到如下结果，P值小于0.05可以认为不同分位数回归的截距之间值是有差异的。
```r
fit1 <- rq(foodexp ~ income, tau = .25, data = engel)
fit2 <- rq(foodexp ~ income, tau = .5, data = engel)
fit3 <- rq(foodexp ~ income, tau = .75, data = engel)
anova(fit1, fit2, fit3)
```
结果有显著差异。

## 混合线性模型lmm
### 固定因子与随机因子

混合线性模型考察既有随机因子，又有固定因子的模型的线性回归问题。而关于固定因子和随机因子，可参考csdn上的一篇博文：[固定效应模型与随机效应模型](https://blog.csdn.net/fjsd155/article/details/94313748)
固定效应和随机效应的选择是大家做面板数据常常要遇到的问题，一个常见的方法是做huasman检验，即先估计一个随机效应，然后做检验，如果拒绝零假设，则可以使用固定效应，反之如果接受零假设，则使用随机效应。但这种方法往往得到事与愿违的结果。另一个想法是在建立模型前根据数据性质确定使用那种模型，比如数据是从总体中抽样得到的，则可以使用随机效应，比如从N个家庭中抽出了M个样本，则由于存在随机抽样，则建议使用随机效应，反之如果数据是总体数据，比如31个省市的Gdp，则不存在随机抽样问题，可以使用固定效应。同时，从估计自由度角度看，由于固定效应模型要估计每个截面的参数，因此随机效应比固定效应有较大的自由度.
***

### 混合线性回归-lme4

混合线性回归：回归是研究因变量与自变量之间的关系，但是如果因变量有着不同的来源呢？比如研究全市学生成绩与学习时间的关系，但是全市众多学校每所学校水平不同，为了平衡/解决学校水平不同的差异，就需要利用混合线性回归。

我们还是使用之前的[加拉帕戈斯地雀模型](http://bioquest.org/birdd/morph.php)作为练习数据，对此进行预处理：
```r {class=line-numbers}
library(tidyverse)
morph <- read.csv("data/raw/Morph_for_Sato.csv", stringsAsFactors = FALSE, strip.white = TRUE)
#使用tolower函数改成小写
names(morph) <- tolower(names(morph))
morph <- morph %>%
  dplyr::select(islandid, taxonorig, genusl69, speciesl69, sex, wingl, beakh, ubeakl) %>%
  dplyr::rename(taxon = taxonorig, genus = genusl69, species = speciesl69)
morph <- data.frame(na.omit(morph))
morph <- dplyr::filter(morph, genus == "Geospiza") %>% as_data_frame()
d <- morph
#把数据保存为RDS格式
saveRDS(d, file = "data/generated/morph-geospiza.rds")
```
我们希望看一看喙长和翅膀长度有什么关系没有：
```r
#取对数缩小范围，并直接在图里进行线性回归，geom_smooth是添加平滑曲线
ggplot(d, aes(log(wingl), log(beakh))) +
  geom_point(aes(colour = taxon)) +
  geom_smooth(method = "lm")
```
![简单的回归](R/Rplot16.jpeg)

但如果我们把`aes(colour = taxon)`移到上面会把我们默认的回归对象变成每个种
我们可以发现在种水平每个种的规律都不一样。
```r
ggplot(d, aes(log(wingl), log(beakh), colour = taxon)) +
  geom_point() +
#se指置信区间
  geom_smooth(method = "lm", se = FALSE)
```
![简单的回归2](R/Rplot17.jpeg)

接下来使用R包“lme4”来对数据进行混合线性回归，我们要建立一个混合效应模型，让每个分类单元都有自己的随机截距：
```r
library(tidyr)
library(Matrix)
library(lme4)
#先普通地线性回归看一看
m_lm <- lm(log(beakh) ~ log(wingl), data = d)
summary(m_lm)
#进行混合线性模型，1表示截距而taxon为随机因子
m_lmer <-  lmer(log(beakh) ~ log(wingl) + (1 | taxon), data = d)
#利用这个混合模型进行预测，再拿预测值画回归线
d$predict_lmer <- predict(m_lmer)
ggplot(d, aes(log(wingl), log(beakh), colour = taxon)) +
  geom_point(alpha = 0.1) + 
  geom_line(aes(y = predict_lmer))
```
我们发现各个元素分别回归的斜率已经相同了。
![混合线性回归](R/Rplot18.jpeg)

我们可以将这些线条合并，从群落结构去考察：
```r
#使用re.form = NA函数来合并taxon数据
d$predict_lmer_population <- predict(m_lmer, re.form = NA)
ggplot(d, aes(log(wingl), log(beakh), colour = taxon)) +
  geom_point(alpha = 0.1) +
  geom_line(aes(y = predict_lmer)) +
  geom_line(aes(y = predict_lmer_population), colour = "black", size = 1)
```
![混合线性回归2](R/Rplot19.jpeg)

我们也可以提取所有这些物种各个的随机截距
```r
ranef(m_lmer)
```
或者查看这些随机截距的均值
```r
round(mean(ranef(m_lmer)$taxon[, 1]), 2)
```
所以所以每个分类单元的截距估计值等于“固定效应”截距加上“随机”偏差。
最后我们可以查看回归效果：
```r
summary(m_lmer)
##关于数据的查看上文已经提及，不再赘述：
#Linear mixed model fit by REML ['lmerMod']
#Formula: log(beakh) ~ log(wingl) + (1 | taxon)
#   Data: d
#
#REML criterion at convergence: -3656.1
#
#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-6.2724 -0.6086 -0.0306  0.6211  3.6835 
#
#Random effects:
# Groups   Name        Variance Std.Dev.
# taxon    (Intercept) 0.050586 0.22491 
# Residual             0.004278 0.06541 
#Number of obs: 1434, groups:  taxon, 15
#
#Fixed effects:
#            Estimate Std. Error t value
#(Intercept) -2.60848    0.18819  -13.86
#log(wingl)   1.18318    0.04232   27.96
#
#Correlation of Fixed Effects:
#           (Intr)
#log(wingl) -0.951
```
一个小技巧，可以通过`arm::display(m_lmer)`(这里使用了arm包）命令来让我们的结果更加直观：
```r
#lmer(formula = log(beakh) ~ log(wingl) + (1 | taxon),data = d)
#            coef.est coef.se
#(Intercept) -2.61     0.19  
#log(wingl)   1.18     0.04  
#
#Error terms:
# Groups   Name        Std.Dev.
# taxon    (Intercept) 0.22    
# Residual             0.07    
#---
#number of obs: 1434, groups: taxon, 15
#AIC = -3648.1, DIC = -3672.8
#deviance = -3664.4 
```
不过需要注意的是，组的数量少于5个时就不建议使用混合线性模型了。

综上，对于lmer()下的混合线性回归，有：
```r
fit = lmer(data = , formula = DV ~ Fixed_Factor + (Random_intercept + Random_Slope | Random_Factor))
```
其中formula为表达式，DV为因变量，Fixed_Factor为固定因子（自变量），Random_intercept为随机截距（不同群体的因变量分布不同，有人在起点，有人在终点），Random_Slope为随机斜率（不同群体受到固定因子的影响不同），Random_Factor为随机因子
。
```r
# nlme::lme()
y ~ x, random = ~1 | group # varying intercepts
y ~ x, random = ~1 + x | group # varying intercepts and slopes
# 上述即为x受到了group的影响（斜率/效果上的）
y ~ x, random = list(group = pdDiag(~ x)) # uncorrelated varying intercepts and slopes
y ~ x, random = ~1 | group/subgroup # nested
# crossed... some structures possible but not easy,
# e.g. Pinheiro and Bates 2000 p163
# e.g. http://stackoverflow.com/questions/36643713/how-to-specify-different-random-effects-in-nlme-vs-lme4

# lme4::lmer()
y ~ x + (1 | group) # varying intercepts
y ~ x + (1 + x | group) # varying intercepts and slopes
y ~ x + (1 | group) + (x - 1 | group) # uncorrelated varying intercepts and slopes
y ~ x + (1 | group/subgroup) # nested
y ~ x + (1 | group1) + (1 | group2) # varying intercepts, crossed
y ~ x + (1 + x | group1) + (1 + x | group2) # varying intercepts and slopes, crossed
```
一群研究生正在尝试将log（青蛙密度）作为植被特征的函数来建模。在以下情况下（使用“lme4::lmer”）它们的模型会是什么样子？只给出具有随机截距的情况，以及具有随机截距和斜率的情况。
Jerry测量了6个独立池塘（“pond”）内1个样带（“transect”）的青蛙密度和植被特征。
```r
# random intercepts
log(frog_dens) ~ vegetation + 
  (1 | pond) # exercise
# random slopes and intercepts:
log(frog_dens) ~ vegetation + 
  (1 + vegetation | pond) # exercise
```
### 混合线性回归-lmerTest
R语言中实现混合线性模型可以使用lme4包或lmerTest包，这里以lmerTest包为例，其基本表达式为：
```r
fit = lmer(data = , formula = DV ~ Fixed_Factor + (Random_intercept + Random_Slope | Random_Factor))
```
其中data为我们要处理的数据集，formula为表达式，DV是因变量，Fixed_Factor是固定因子（自变量），Random_intercept是随机截距（可以理解为因变量分布的不同？），Random_Slope是随机斜率，即认为不同群体受固定因子的影响不同，Random_Factor是随机因子。
我们以politeness数据为例进行计算：
politeness数据可以在github：<https://github.com/usplos/Eye-movement-related/blob/master/politeness_data.csv>中获得，本篇关于混合线性的模型的计算也源自该项目。该数据收集了若干被试（subject）的性别（gender），以及用不同的态度（attitude）在不同场合（scenario）下说话的音高（frequency）。 这是一个典型的被试内设计（7 * 2设计）。
先打开数据并加载相关r包：
```r
politeness = readr::read_csv('/Users/desktop/r/politeness_data.csv')
library(lme4)
library(Matrix)
library(lmerTest)
#将scenairo变为因子型变量（离散型，原来是字符型）
politeness$scenario = factor(politeness$scenario)
politeness
```
进行混合线性计算：
```r
#音高与固定因子态度和场合的关系，随机因子是性别与被试，它们基于的设计矩阵是态度，即为我们认为性别和被试对回归的影响的随机的
fit1 = lmer(frequency ~ scenario * attitude + attitude|subject + attitude|gender, data = politeness)
```
这时候系统提示：`boundary (singular) fit: see ?isSingular`这表明有些效应是彼此的线性组合或者某个地方的方差是0。当p>>n(比样本更多的参数)时也会发生这种情况。这也意味着模型可能被过度拟合和/或遭受数值稳定性问题。
```r
summary(fit1)
##得到结果：
#Linear mixed model fit by REML. t-tests use #Satterthwaite's method ['lmerModLmerTest']
#Formula: frequency ~ scenario * attitude + (1 + attitude | subject) +  
#    (1 + attitude | gender)
#   Data: politeness
#
#REML criterion at convergence: 680.1
#
#Scaled residuals: 
#     Min       1Q   Median       3Q      Max 
#-1.65342 -0.68640 -0.03677  0.50256  2.85422 
##这里是随机效应的结果，variance为方差，Std.Dev为标准差，可以看出确实对不同的被试组合性别组而言，态度的影响是不同的
#Random effects:
# Groups   Name        Variance  Std.Dev. Corr 
# subject  (Intercept) 6.037e+02 24.5696       
#          attitudepol 1.076e-02  0.1037  1.00 
# gender   (Intercept) 6.467e+03 80.4167       
#          attitudepol 1.118e+02 10.5749  -1.00
# Residual             6.101e+02 24.7001       
#Number of obs: 83, groups:  subject, 6; gender, 2
##固定效应的结果如下，我们发现场合3和4的音高是较显著的
#Fixed effects:
#                      Estimate Std. Error      df t value Pr(>|t|)   
#(Intercept)            180.767     58.615   1.065   3.084  0.18720   
#scenario2               17.450     14.261  63.998   1.224  0.22557   
#scenario3               46.667     14.261  63.998   3.272  0.00172 **
#scenario4               44.833     14.261  63.998   3.144  0.00253 **
#scenario5               16.800     14.261  63.998   1.178  0.24313   
#scenario6                8.867     14.261  63.998   0.622  0.53631   
#scenario7               18.133     14.261  63.998   1.272  0.20813   
#attitudepol             -9.717     16.102   9.583  -0.603  0.56023   
#scenario2:attitudepol   15.133     20.168  63.998   0.750  0.45578   
#scenario3:attitudepol  -31.283     20.168  63.998  -1.551  0.12579   
#scenario4:attitudepol   -4.650     20.168  63.998  -0.231  0.81839   
#scenario5:attitudepol   -4.783     20.168  63.998  -0.237  0.81327   
#scenario6:attitudepol  -14.703     20.701  64.030  -0.710  0.48011   
#scenario7:attitudepol  -30.033     20.168  63.998  -1.489  0.14135   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Correlation matrix not shown by default, as p = 14 > 12.
#Use print(x, correlation=TRUE)  or
#    vcov(x)        if you need it
#
#convergence code: 0
#boundary (singular) fit: see ?isSingular
```
但是这里的固定效应不是主效应和交互作用，要查看主效应和交互作用需要用anova()函数得到；
```r
anova(fit1)
#Type III Analysis of Variance Table with #Satterthwaite's method
#                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
#scenario          19400.1  3233.4     6 64.006  5.2998 #0.000173 ***
#attitude           2789.7  2789.7     1  1.143  4.5725 0.253068    
#scenario:attitude  4985.4   830.9     6 64.006  1.3619 0.243577    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
可见只有场景的主效应显著，态度的主效应和交互作用都不显著。在上面建立的模型中，包含随机斜率和随机截距，但是注意到，有两个固定效应，是把两个固定效应及其交互作用全都作为随机效应，还是选其中部分作为随机效应呢？这里我们课题组的标准是：首先考虑全模型，即如下命令：
```r
fitAll = lmer(frequency ~ scenario * attitude + (attitude * scenario|subject) + (attitude * scenario|gender), data = politeness)
```
然后报错，观测量的个数少于随机因子的个数，因此移除交互作用：
```r
fitAll2 = lmer(frequency ~ scenario * attitude + (attitude + scenario|subject) + (attitude + scenario|gender), data = politeness)
```
结果模型似乎炸了，变成金字塔了，于是我们尝试移除一些随机因子：同时设立一个只有随机截距的零模型：
```r
fitAll3_1 = lmer(frequency ~ scenario * attitude + (attitude|subject) + (attitude + scenario|gender), data = politeness);
fitAll3_2 = lmer(frequency ~ scenario * attitude + (scenario|subject) + (attitude + scenario|gender), data = politeness);
fitAll3_3 = lmer(frequency ~ scenario * attitude + (attitude+ scenario|subject) + (attitude|gender), data = politeness);
fitAll3_4 = lmer(frequency ~ scenario * attitude + (attitude + scenario|subject) + (scenario|gender), data = politeness)
fitZero = lmer(frequency ~ scenario * attitude + (1|subject) + (1|gender), data = politeness)
```
使用`anova()`分别比较各个模型和零模型的P值，发现都不显著，这时选取P值最小的作为实际在论文中使用的模型，即选取fitAll3_3。

### 实战：拟合、解释并检验lmm

这里还是使用gapminder数据集
读取数据：
```r
library(lme4)
library(tidyverse)
theme_set(theme_light())
gap <- gapminder::gapminder
gap <- mutate(gap, decade = (year - 1990)/10,
  log_gdp_percap_cent = log(gdpPercap) - mean(log(gdpPercap)))
```
我们增加了以1990年为中心的十年新列，并增加了一个“log_gdp_percap_cent”列，即以平均值为中心的人均gdp的对数列，这样可以提升结果的准确性。
在开始建模之前，让我们先画一个数据图。
```r
ggplot(gap, aes(log_gdp_percap_cent, lifeExp, colour = continent)) + 
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = country), method = "lm", se = FALSE, alpha = 0.5, lwd = 0.5)
```
![优势比](R/Rplot125.jpeg)
我们已经将一个更简单的模型拟合到这个数据集。让我们试着拟合一个线性混合效应模型（使用`lmer()`），根据“log_gdp_percap_cent”和“decade”预测“lifeExp”。根据上面的情节和你对世界的理解，决定是只允许截距变化还是截距和斜率因国家而异。
```r
m <- lmer(lifeExp ~
    log_gdp_percap_cent + decade + # exercise
    (log_gdp_percap_cent + decade | country), data = gap) # exercise
arm::display(m)
# lmer(formula = lifeExp ~ log_gdp_percap_cent + decade + (log_gdp_percap_cent + 
#     decade | country), data = gap)
#                     coef.est coef.se
# (Intercept)         62.02     0.64  
# log_gdp_percap_cent  4.29     0.38  
# decade               2.56     0.18  

# Error terms:
#  Groups   Name                Std.Dev. Corr        
#  country  (Intercept)         7.09                 
#           log_gdp_percap_cent 3.21     -0.17       
#           decade              2.01      0.01 -0.72 
#  Residual                     1.96                 
# ---
# number of obs: 1704, groups: country, 142
# AIC = 8410.1, DIC = 8387.2
# deviance = 8388.7
```
国家在各大洲之间相互依存。尝试拟合上述模型的另一个版本，其中随机效果因大陆和嵌套在大陆中的国家而异：
```r
m2 <- lmer(lifeExp ~ 
    log_gdp_percap_cent + decade + # exercise
    (log_gdp_percap_cent + decade | continent/country), data = gap) # exercise
arm::display(m2)
# lmer(formula = lifeExp ~ log_gdp_percap_cent + decade + (log_gdp_percap_cent + 
#     decade | continent/country), data = gap)
#                     coef.est coef.se
# (Intercept)         64.92     2.84  
# log_gdp_percap_cent  3.38     0.46  
# decade               2.59     0.43  

# Error terms:
#  Groups            Name                Std.Dev.
#  country:continent (Intercept)         5.14    
#                    log_gdp_percap_cent 2.66    
#                    decade              1.66    
#  continent         (Intercept)         5.97    
#                    log_gdp_percap_cent 0.59    
#                    decade              0.86    
#  Residual                              1.97    
#  Corr        
             
#  -0.14       
#   0.04 -0.67 
             
#  -1.00       
#  -0.36  0.36 
             
# ---
# number of obs: 1704, groups: country:continent, 142; continent, 5
# AIC = 8317.9, DIC = 8292.6
# deviance = 8289.2
```
看看你拟合的第一个模型的拟合值和残差。
你可能想把它们按大陆分开，让它们更容易看。可以使用“broom:：augment()`并使用ggplot自己绘制它们，也可以像我们在方差结构练习中使用的那样使用快捷语法。
还可以根据模型中的预测值和数据集中未包含在模型中的任何其他投影仪绘制残差。你觉得这些怎么样？
```r
plot(m, resid(.) ~ fitted(.) | continent, abline = 0)
plot(m, resid(.) ~ log(pop) | continent, abline = 0) # exercise
plot(m, resid(.) ~ decade | continent, abline = 0) # exercise
plot(m, resid(.) ~ log_gdp_percap_cent | continent, abline = 0) # exercise
```
![拟合](R/Rplot126.jpeg)

检查随机效应是否近似正态分布：
```r
re <- ranef(m)
qqnorm(re$country[,"(Intercept)"])
qqline(re$country[,"(Intercept)"])

qqnorm(re$country[,"log_gdp_percap_cent"])
qqline(re$country[,"log_gdp_percap_cent"])
```
尝试绘制模型预测。让我们关注人均GDP的影响。因此，我们必须将预测值“decade”设置为某个值。我们把它设为'0'，它代表1990年，因为我们是通过减去1990来计算这个列的。
```r
newdata <- mutate(gap, decade = 0)
newdata$predict <- predict(m, newdata = newdata)
# 绘图
ggplot(newdata, aes(log_gdp_percap_cent, lifeExp, colour = continent)) + 
  geom_line(aes(y = predict, group = country), alpha = 0.5)
```
![预测](R/Rplot127.jpeg)

让我们添加一个组级预测器，代表每个国家人均GDP的平均对数。同时，让我们把每个国家人均GDP中心化：
```r
gap <- group_by(gap, country) %>% 
  mutate(mean_log_gdp_percap_cent = mean(log_gdp_percap_cent), 
    log_gdp_percap_group_cent = log_gdp_percap_cent - mean(log_gdp_percap_cent)) %>%
  ungroup()
```
现在将其添加到您适合的初始模型中。这个模型不会收敛于大陆内部的嵌套随机效应。所以，让影响因国家而异:
```r
m3 <- lmer(lifeExp ~ log_gdp_percap_group_cent + decade + 
    mean_log_gdp_percap_cent + 
    (log_gdp_percap_cent + decade |  country), data = gap) 
# 绘图
newdata <- mutate(gap, decade = 0)
newdata$predict3 <- predict(m3, newdata = newdata)
ggplot(newdata, aes(log_gdp_percap_cent, lifeExp, colour = continent)) + 
  geom_line(aes(y = predict3, group = country), alpha = 0.5)
```
![预测](R/Rplot128.jpeg)

和上一个模型对比一下：
```r
arm::display(m)
# lmer(formula = lifeExp ~ log_gdp_percap_cent + decade + (log_gdp_percap_cent + 
#     decade | country), data = gap)
#                     coef.est coef.se
# (Intercept)         62.02     0.64  
# log_gdp_percap_cent  4.29     0.38  
# decade               2.56     0.18  

# Error terms:
#  Groups   Name                Std.Dev. Corr        
#  country  (Intercept)         7.09                 
#           log_gdp_percap_cent 3.21     -0.17       
#           decade              2.01      0.01 -0.72 
#  Residual                     1.96                 
# ---
# number of obs: 1704, groups: country, 142
# AIC = 8410.1, DIC = 8387.2
# deviance = 8388.7
arm::display(m3)
# lmer(formula = lifeExp ~ log_gdp_percap_group_cent + decade + 
#     mean_log_gdp_percap_cent + (log_gdp_percap_cent + decade | 
#     country), data = gap)
#                           coef.est coef.se
# (Intercept)               62.18     0.54  
# log_gdp_percap_group_cent  3.05     0.40  
# decade                     2.77     0.18  
# mean_log_gdp_percap_cent   8.15     0.51  

# Error terms:
#  Groups   Name                Std.Dev. Corr        
#  country  (Intercept)         5.86                 
#           log_gdp_percap_cent 3.05     -0.23       
#           decade              1.93      0.20 -0.68 
#  Residual                     1.95                 
# ---
# number of obs: 1704, groups: country, 142
# AIC = 8341.7, DIC = 8317.1
# deviance = 8318.4
```
接下来进行模型解释
```r
round(fixef(m3)["decade"], 1)
# decade 
#    2.8
# 查看每十年预期寿命增长的95%置信区间 
fe <- fixef(m3)["decade"] # exercise
se <- arm::se.fixef(m3)["decade"] # exercise
fe + c(-1.96, 1.96) * se # exercise
# [1] 2.421557 3.110131
```
人均国内生产总值的影响是在国家内部还是在国家之间哪个更强？
```r
arm::display(m3)
# lmer(formula = lifeExp ~ log_gdp_percap_group_cent + decade + 
#     mean_log_gdp_percap_cent + (log_gdp_percap_cent + decade | 
#     country), data = gap)
#                           coef.est coef.se
# (Intercept)               62.18     0.54  
# log_gdp_percap_group_cent  3.05     0.40  
# decade                     2.77     0.18  
# mean_log_gdp_percap_cent   8.15     0.51  

# Error terms:
#  Groups   Name                Std.Dev. Corr        
#  country  (Intercept)         5.86                 
#           log_gdp_percap_cent 3.05     -0.23       
#           decade              1.93      0.20 -0.68 
#  Residual                     1.95                 
# ---
# number of obs: 1704, groups: country, 142
# AIC = 8341.7, DIC = 8317.1
# deviance = 8318.4
```
除去国内生产总值的差异，加拿大、中国和津巴布韦每十年预期寿命的估计变化是多少?(提示:要么使用' coef() '，要么结合使用' fixef() '和' ranef() '的固定和随机效果。看看输出，这里不需要用R代码提取值。)
```r
re <- coef(m3)$country # exercise
re$country <- row.names(re) # exercise
filter(re, country == "Canada") %>% pull(decade) 
# [1] 1.497968
filter(re, country == "China") %>% pull(decade) 
# [1] 6.883605
filter(re, country == "Zimbabwe") %>% pull(decade) 
# [1] -1.362036
```


## 地理自相关——基于nlme

学习通过nlme中的相关结构识别和处理空间自相关。我们将使用nlme包中的一个示例数据集。这个数据集代表了在不同纬度和经度采集的不同小麦品种的产量。
```r
library(nlme)
library(tidyverse)
d <- nlme::Wheat2
glimpse(d)
ggplot(d, aes(longitude, latitude, size = yield, colour = variety)) + 
  geom_point()
```
![地理分布状况](R/Rplot96.jpeg)

到目前为止，我们已经使用nlme包中的lme函数来拟合线性混合效应模型。这个软件包还具有函数“gls”，它允许您在拟合线性模型时使用nlme的特性而不产生随机效果。
```r
m1 <- gls(yield ~ variety - 1, data = d)
m1
# Generalized least squares fit by REML
#   Model: yield ~ variety - 1 
#   Data: d 
#   Log-restricted-likelihood: -620.3709

# Coefficients:
#   varietyARAPAHOE      varietyBRULE 
#           29.4375           26.0750 
#   varietyBUCKSKIN    varietyCENTURA 
#           25.5625           21.6500 
#  varietyCENTURK78   varietyCHEYENNE 
#           30.3000           28.0625 
#       varietyCODY       varietyCOLT 
#           21.2125           27.0000 
#       varietyGAGE  varietyHOMESTEAD 
#           24.5125           27.6375 
#   varietyKS831374     varietyLANCER 
#           24.1250           28.5625 
#    varietyLANCOTA    varietyNE83404 
#           26.5500           27.3875 
#    varietyNE83406    varietyNE83407 
#           24.2750           22.6875 
#    varietyNE83432    varietyNE83498 
#           19.7250           30.1250 
#    varietyNE83T12    varietyNE84557 
#           21.5625           20.5250 
#    varietyNE85556    varietyNE85623 
#           26.3875           21.7250 
#    varietyNE86482    varietyNE86501 
#           24.2875           30.9375 
#    varietyNE86503    varietyNE86507 
#           32.6500           23.7875 
#    varietyNE86509    varietyNE86527 
#           26.8500           22.0125 
#    varietyNE86582    varietyNE86606 
#           24.5375           29.7625 
#    varietyNE86607   varietyNE86T666 
#           29.3250           21.5375 
#    varietyNE87403    varietyNE87408 
#           25.1250           26.3000 
#    varietyNE87409    varietyNE87446 
#           21.3750           27.6750 
#    varietyNE87451    varietyNE87457 
#           24.6125           23.9125 
#    varietyNE87463    varietyNE87499 
#           25.9125           20.4125 
#    varietyNE87512    varietyNE87513 
#           23.2500           26.8125 
#    varietyNE87522    varietyNE87612 
#           25.0000           21.8000 
#    varietyNE87613    varietyNE87615 
#           29.4000           25.6875 
#    varietyNE87619    varietyNE87627 
#           31.2625           23.2250 
#     varietyNORKAN    varietyREDLAND 
#           24.4125           30.5000 
# varietyROUGHRIDER    varietySCOUT66 
#           21.1875           27.5250 
#  varietySIOUXLAND     varietyTAM107 
#           30.1125           28.4000 
#     varietyTAM200       varietyVONA 
#           21.2375           23.6000 

# Degrees of freedom: 224 total; 168 residual
# Residual standard error: 7.711373
```
我们可以自己提取和绘制空间上的残差（通过颜色）。
```r
d$res <- as.numeric(residuals(m1))
ggplot(d, aes(longitude, latitude, colour = res)) + 
  geom_point(size = 5) + scale_color_gradient2()
```
![地理分布残差](R/Rplot97.jpeg)

我们可以清楚地看到，残差具有空间聚集模式。为什么这是个问题？
nlme软件包附带了一个用于绘制半变异函数的内置函数。这表示在彼此距离增加的情况下，残差之间的平均平方差的一半。所以这和相关性成反比。小值意味着残差彼此非常相似。这是空间统计中的一个常见情节，你可以在网上找到很多关于这个主题的参考资料。
```r
plot(Variogram(m1, form = ~ latitude + longitude, data = d))
```
![残差相关](R/Rplot98.jpeg)

如果我们看半变异函数，我们可以猜测一个范围值和掘金效应的良好初始值。对于“nlme:：corsphere()`相关结构，nugget表示截距值，range表示半变异函数达到1的距离。
通过观察上面的半变异函数，我们可以看到这些值应该是什么。我们将给相关函数的起始值接近我们所需要的
期待。
```r
m2 <- update(m1, 
  corr = corSpher(c(30, 0.2), form = ~ latitude + longitude, nugget = TRUE))
m2
```
事实上，这是一个例子，如果我们不给相关结构适当的起始值，它会得到错误的答案：
```r
m3 <- update(m1, 
  corr = corSpher(form = ~ latitude + longitude, nugget = TRUE))
m3
```
让我们试着在空间上绘制残差，并制作另一个半变异函数。请再次注意，为了将相关结构合并到残差计算中，使用`type=“normalized”`非常重要。
```r
d$res2 <- as.numeric(residuals(m2, type = "normalized"))
ggplot(d, aes(longitude, latitude, colour = res2)) + geom_point(size = 5) +
  scale_color_gradient2()
plot(Variogram(m2, form = ~ latitude + longitude, resType = "normalized", data = d))
```
![空间自相关](R/Rplot99.jpeg)

![空间自相关](R/Rplot100.jpeg)

我们还可以查看他们的deltaAIC：
```r
bbmle::AICtab(m1, m2)
#    dAIC  df
# m2   0.0 59
# m1 169.3 57
```
这种方法建模空间相关性的一个缺点是，没有简单的方法来提取预测的空间曲面。另一种处理空间自相关的方法是用GAM将空间过程建模为二维光滑项。GAM超出了本次研讨会的范围，但重要的是要知道它们的存在。
GAM和GLMs一样，只是可以允许预测因子沿着平滑的摆动线。在拟合算法中可以客观地确定摆动程度。
作为演示，我们将用GAM安装此模型的一个版本。在R中安装这些模型的主要包是mgcv。我们将用最大似然法修正“gls”模型，以便与GAM进行比较。
```r
library(mgcv)
m1_ml <- gls(yield ~ variety - 1, data = d, method = "ML")
m_gam1 <- gam(yield ~ variety - 1, data = d) # the same
bbmle::AICtab(m_gam1, m1_ml)
m_gam2 <- gam(yield ~ variety - 1 + te(latitude, longitude), data = d)

bbmle::AICtab(m_gam1, m_gam2)
#        dAIC df
# m_gam1  0   57
# m1_ml   0   57
#        dAIC  df  
# m_gam2   0.0 75.2
# m_gam1 269.8 57  
```
我们现在可以绘制出一个空间预测图，并像以前一样在空间上检查残差：：
```r
plot(m_gam2, pers = TRUE)
```
![空间预测](R/Rplot101.jpeg)

![空间残差](R/Rplot102.jpeg)

### 数据的分布与回归
回顾GLMs和GLMMs中最常用的概率分布
了解参数如何影响分布形状
正态分布
```r
library(manipulate)

x <- seq(-10, 10, length.out = 300)
manipulate({
  y <- dnorm(
    x = x,
    mean = mean,
    sd = sd)
  plot(x = x, y = y, type = "l", ylim = c(0, max(y)),
    main = "Normal")},
  mean = slider(-10, 10, 0),
  sd = slider(0.1, 10, 2))
```
![分布](R/Rplot103.jpeg)

log变换下的正态分布
```r
x <- seq(0, 100, length.out = 300)
manipulate({
  y <- dlnorm(
    x = x,
    meanlog = meanlog,
    sdlog = sdlog)
  plot(x = x, y = y, type = "l", xlim = c(0, 100), ylim = c(0, max(y)),
    main = "Lognormal")},
  meanlog = slider(0.1, 10, 3),
  sdlog = slider(0.01, 2, 1))
```
![分布](R/Rplot104.jpeg)

Gamma分布
可以理解为n个指数分布的独立随机变量的加总
```r
x <- seq(0, 100, length.out = 300)
manipulate({
  y <- dgamma(
    x = x,
    shape = shape,
    scale = scale)
  plot(x = x, y = y, type = "l", ylim = c(0, max(y)),
    main = "Gamma")},
  shape = slider(1, 10, 3),
  scale = slider(0.1, 30, 2))
```
![分布](R/Rplot105.jpeg)

接下来是离散分布
泊松分布
```r
x <- seq(0, 100)
manipulate({
  y <- dpois(x,
    lambda = lambda)
  plot(x = x, y = y, type = "h", xlim = c(0, 60), ylim = c(0, max(y)),
    main = "Poisson")},
  lambda = slider(0.1, 30, 3))
```
![分布](R/Rplot106.jpeg)

负二项分布
```r
x <- seq(0, 100)
manipulate({
  y <- dnbinom(
    x = x,
    size = size,
    mu = mu)
  plot(x = x, y = y, type = "h", ylim = c(0, max(y)),
    main = "Negative binomial")},
  size = slider(0.1, 10, 1),
  mu = slider(0.1, 60, 4))
```
![分布](R/Rplot107.jpeg)

二项分布
```r
x <- seq(0, 10)
manipulate({
  y <- dbinom(x,
    size = size,
    prob = prob)
  plot(x = x, y = y, type = "h", xlim = c(0, 10), ylim = c(0, max(y)),
    main = "Binomial")},
  size = slider(1, 10, initial = 1, step = 1),
  prob = slider(0, 1, initial = 0.5, step = 0.1))
```
size = 1, prob = 0.5
![分布](R/Rplot108.jpeg)

size = 10, prob = 0.4
![分布](R/Rplot109.jpeg)

## 广义线性模型glm

在R语言中可以通过`glm()`函数解决广义线性模型，此处我们运用logistic模型进行广义线性回归：
```r
glm(formula,family = gaussian,data,...)
```
其中formula为要拟合的模型，而family为分布族，包括正态分布、泊松分布、二项分布等...分布族可以通过link=来选择连接函数，data为数据框。
首先导入数据，数据来源于[我的github](https://github.com/Vendredii/Rstats)，最初源于教材《多元线性回归与R》：
```r
data = readr::read_csv('/Users/desktop/r/result.csv')
#其中x1为视力状况，1是好0是不好；x2为年龄；x3为驾车教育，1是有0是没有；y为是否出过事故，1是有0是没有
data
```
在这里y是因变量，只有两个值，因此我们可以把它看作是成功概率为$p$的Bernoulli试验的结果（这种方法在进行二分类问题时很有用！），现在用Logistic回归模型进行分析：
假定模型为：
$$\textup{ln}(\frac{p}{1-p})=\beta _0+\beta _1x_1+\beta _2x_2+\beta _3x_3$$
有：
```r
logit.glm <- glm(y ~ X1 + X2 + X3, family = binomial,data = data)
summary(logit.glm)
#glm(formula = y ~ X1 + X2 + X3, family = binomial, data = data)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.5636  -0.9131  -0.7892   0.9637   1.6000  
#
#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)  
#(Intercept)  0.597610   0.894831   0.668   0.5042  
#X1          -1.496084   0.704861  -2.123   0.0338 *
#X2          -0.001595   0.016758  -0.095   0.9242  
#X3           0.315865   0.701093   0.451   0.6523  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 62.183  on 44  degrees of freedom
#Residual deviance: 57.026  on 41  degrees of freedom
#AIC: 65.026
#
#Number of Fisher Scoring iterations: 4
```
由此得到初步的Logistic回归模型：
$$p=\frac{\textup{exp}(0.598-1.496x_1-0.002x_2+0.316x_3)}{1+\textup{exp}(0.598-1.496x_1-0.002x_2+0.316x_3)}$$
即为：
$$\textup{Logit}(p)=0.598-1.496x_1-0.002x_2+0.316x_3$$
由于参数$\beta _2$和$\beta _3$没有通过P值检验，可通过`step()`作变量筛选，不断筛选出AIC值最小的结果：
```r
logit.step <- step(logit.glm, direction = "both")
summary(logit.step)
#glm(formula = y ~ X1, family = binomial, data = data)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.4490  -0.8782  -0.8782   0.9282   1.5096  
#
#Coefficients:
#            Estimate Std. Error z value Pr(>|z|)  
#(Intercept)   0.6190     0.4688   1.320   0.1867  
#X1           -1.3728     0.6353  -2.161   0.0307 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 62.183  on 44  degrees of freedom
#Residual deviance: 57.241  on 43  degrees of freedom
#AIC: 61.241
#
#Number of Fisher Scoring iterations: 4
```
可以得到新的回归方程为：
$$p=\frac{\textup{exp}(0.619-1.373x_1)}{1+\textup{exp}(0.619-1.373x_1)}$$
可见$p_1=0.32,p_2=0.65$，即实力有问题的司机发生交通事故的概率是正常司机的两倍！

### 解决泊松与二项分布中的过度分散

在拟合中可能出错的一个常见问题是，数据中的可变性比分布所允许的要大。对于正态分布、Gamma分布或负二项分布来说，这不是问题，因为这些分布有一个参数，使它们可以根据需要窄或宽。
但是一些分布，特别是泊松分布和二项式分布（有多个样本），对于给定的平均值假设了一个固定的可变性水平。但现实世界很混乱，情况并非总是如此。让我们看看这意味着什么。
我们将从生成已知泊松分布过大的计数数据开始。
```r
library(ggplot2)
library(dplyr)
set.seed(111)
N <- 500
x <- runif(N, -1, 1)
a <- 0.5
b <- 1.3
d <- data_frame(x = x)
inverse_logit <- function(x) plogis(x)

y_true <- exp(a + b * x)

rqpois <- function (n, lambda, d = 1) { # generate random quasipoisson values
  if (d == 1)
    rpois(n, lambda)
  else
    rnbinom(n, size = (lambda / (d - 1)), mu = lambda)
}

set.seed(1234)
y <- rqpois(N, lambda = y_true, d = 5)
plot(x, y)
```
![分布](R/Rplot110.jpeg)

让我们看看我们刚刚创建的数据。
在下面的图中，虚线表示一对一的线（泊松），蓝线表示方差与平均值成线性比例，而不是一对一（此处的真实关系为准泊松），红线表示方差与平均值成二次比例（负二项式）。
为了绘制这个图，我把x轴上的值分成了15个箱子。
```r
d$y <- y
d$x_group <- findInterval(d$x, seq(min(d$x), max(d$x), length.out = 15))
 group_by(d, x_group) %>%
  summarise(m = mean(y), v = var(y)) %>%
  ggplot(aes(m, v)) +
  geom_smooth(method = "lm", 
    formula = y ~ x - 1, se = F, colour = "blue") +
  geom_smooth(method = "lm", 
    formula = y ~ I(x^2) + offset(x) - 1, colour = "red", se = F) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_point()
```
![拟合](R/Rplot111.jpeg)

让我们用泊松分布和日志链接拟合GLM，即使我们知道数据被过度显示。
```r
(m_poisson <- glm(y ~ x, family = poisson(link = "log"), data = d))
# Call:  glm(formula = y ~ x, family = poisson(link = "log"), data = d)

# Coefficients:
# (Intercept)            x  
#      0.4951       1.2743  

# Degrees of Freedom: 499 Total (i.e. Null);  498 Residual
# Null Deviance:	    2274 
# Residual Deviance: 1820 	AIC: 2566
ggeffects::ggpredict(m_poisson, "x", full.data = TRUE) %>%
  ggplot(aes(x, predicted)) +
  geom_line(colour = "red") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point()
```
![拟合](R/Rplot112.jpeg)

如果我们看残差，它们应该和预测的平均值保持不变。在残差中很难看到这些模式。虽然我们知道有过度分散，但这里没有太多可看的东西。
```r
plot(fitted(m_poisson), residuals(m_poisson))
```
![残差](R/Rplot113.jpeg)

可以使用离散验证：
```r
AER::dispersiontest(m_poisson)
# 	Overdispersion test

# data:  m_poisson
# z = 6.0335, p-value = 8.023e-10
# alternative hypothesis: true dispersion is greater than 1
# sample estimates:
# dispersion 
#   4.805502
```
p < 0.01，显著离散。
为了解决过度分散的问题，我们将用拟泊松族对模型进行修正。这只是估计数据的过度分散程度，并适当地调整参数估计的标准误差。
但这意味着我们离开了可能性的世界，不能简单地计算AIC之类的值。（有qAIC等）
```r
(m_qp <- glm(y ~ x, family = quasipoisson(link = "log"), data = d))
# Call:  glm(formula = y ~ x, family = quasipoisson(link = "log"), data = d)

# Coefficients:
# (Intercept)            x  
#      0.4951       1.2743  

# Degrees of Freedom: 499 Total (i.e. Null);  498 Residual
# Null Deviance:	    2274 
# Residual Deviance: 1820 	AIC: NA
confint(m_qp)
#                 2.5 %    97.5 %
# (Intercept) 0.3219956 0.6566294
# x           1.0031108 1.5555660
confint(m_poisson)
#                 2.5 %    97.5 %
# (Intercept) 0.4176408 0.5701738
# x           1.1493870 1.4012934
```
如果我们有重复的试验，我们可以得到二项分布的过度分散的数据。（如果我们是模拟单次试验，即0和1，就不存在过度分散的问题。）
什么时候会这样？举个例子，也许你是在测量青蛙在一个给定的水箱中存活的比例。
在这个例子中，假设每个水箱有30只青蛙，共40个水箱。
让我们模拟一下，在数据过度分散的情况下，以及在没有过度分散的情况下，一些实验后存活下来的青蛙的比例：
```r
set.seed(1)
n <- 30
y <- emdbook::rbetabinom(40, 0.5, size = n, theta=1)
y2 <- rbinom(40, 0.5, size = n)
par(mfrow = c(2, 1))
plot(table(y/n)/length(y), xlim = c(0, 1), ylab = "prop.", 
  main = "Overdispersed")
plot(table(y2/n)/length(y2), xlim = c(0, 1), ylab = "prop.",
  main = "Not overdispersed")
```
![分布](R/Rplot114.jpeg)

![分布](R/Rplot115.jpeg)

我们在这里看到的是每个水箱中存活的青蛙比例的柱状图。请注意，与纯二项分布相比，过度分散场景中的值要分散得多。
让我们用二项误差分布拟合的GLM和准多项式分布允许过度分布的GLM来绘制估计的平均存活比例。
```r
par(mfrow = c(1, 1))
plot(table(y/n)/length(y), xlim = c(0, 1), ylab = "prop.", col = "grey80")
abline(v = 0.5, col = "black", lwd = 10)

ss <- rep(n, length(y))
m <- glm(y/n ~ 1, family = binomial(link = "logit"),
  weights = ss)
ci <- inverse_logit(confint(m))
abline(v = ci, col = "red", lwd = 5)

m2 <- glm(y/n ~ 1, family = quasibinomial(link = "logit"),
  weights = rep(n, length(y)))
ci2 <- inverse_logit(confint(m2))

abline(v = ci2, col = "blue", lwd = 5)
```
![分布](R/Rplot116.jpeg)

在上面的图中，真实值由粗黑色垂直线表示。
二项式GLM 95%置信区间用红色垂直线表示。
拟多项式GLM 95%置信区间用蓝色垂直线表示。
注意，如果不考虑过度分散，我们的置信区间看起来太小了。
由于这是一个关于GLMMs的课程，处理过度分散的另一种方法是为每个水箱建立一个随机拦截模型。但我们现在还不涉及这个。

### 泰坦尼克号乘客的幸存预测
以下数据集来自[Kaggle](https://www.kaggle.com/c/titanic/data)。它代表了泰坦尼克号上的乘客，无论他们是否幸存，以及他们的一些特征。我们将使用代表乘客年龄、乘客性别（女性=1，男性=0）以及他们为机票支付的票价的列。我们将使用这些特征来预测乘客是否幸存（是=1，否=0）。
```r
library(tidyverse)
d <- read_csv("data/raw/titanic.csv")
d <- mutate(d, female = ifelse(Sex == "female", 1, 0))
names(d) <- tolower(names(d))
d <- select(d, survived, age, fare, female) %>% na.omit %>% as_data_frame()
d
# View(d)
head(d)
# A tibble: 6 x 4
#   survived   age  fare female
#      <dbl> <dbl> <dbl>  <dbl>
# 1        0    22  7.25      0
# 2        1    38 71.3       1
# 3        1    26  7.92      1
# 4        1    35 53.1       1
# 5        0    35  8.05      0
# 6        0    54 51.9       0
```
先花几分钟的时间以图形方式浏览数据。当你对这些数据建模时，你看到了什么样的模式，你会期望什么？
年龄与幸存情况
```r
ggplot(d, aes(age, survived, colour = as.factor(female), size = fare)) + # exercise
  geom_point(position = position_jitter(height = 0.2)) # exercise
```
![分布](R/Rplot117.jpeg)

票价与幸存情况
```r
ggplot(d, aes(fare, survived, colour = as.factor(female), size = age)) +  # exercise
  geom_point(position = position_jitter(height = 0.2)) # exercise
```
![分布](R/Rplot118.jpeg)

性别与幸存情况
```r
ggplot(d, aes(age, survived, colour = log(fare))) +  # exercise
  geom_point(position = position_jitter(height = 0.2)) + # exercise
  facet_wrap(~female) # exercise
```
![分布](R/Rplot119.jpeg)

从一个简单的模型开始，有3个预测因子，没有交互作用。我们正在处理二进制数据作为响应。什么样的分布和联系才有意义？
```r
m <- glm(survived ~ 
    age + fare + female, data = d, family = binomial(link = "logit")) # exercise
arm::display(m)
# glm(formula = survived ~ age + fare + female, family = binomial(link = "logit"), 
#     data = d)
#             coef.est coef.se
# (Intercept) -1.41     0.23  
# age         -0.01     0.01  
# fare         0.01     0.00  
# female       2.35     0.19  
# ---
#   n = 714, k = 4
#   residual deviance = 716.1, null deviance = 964.5 (difference = 248.4)
plot(ggeffects::ggpredict(m)) %>%
  cowplot::plot_grid(plotlist = .)

sjPlot::plot_model(m, type = "est")
```
![预测曲线](R/Rplot120.jpeg)

![优势比](R/Rplot121.jpeg)

请注意，`sjPlot::plot_model()`显示的是*优势比*，而不是我们在`summary()`或`arm::display()``中看到的*对数*优势比。
在这个初始模型中，女性乘客的存活几率比男性乘客大多少？
现在尝试添加所有双向交互。记住这是有捷径的。尝试通过AIC将此模型与没有任何交互的模型进行比较。
```r
m2 <- glm(survived ~
    (age + fare + female) ^ 2, data = d, family = binomial(link = "logit")) # exercise
arm::display(m2)
# glm(formula = survived ~ (age + fare + female)^2, family = binomial(link = "logit"), 
#     data = d)
#             coef.est coef.se
# (Intercept) -0.66     0.35  
# age         -0.03     0.01  
# fare         0.00     0.01  
# female       0.98     0.45  
# age:fare     0.00     0.00  
# age:female   0.04     0.01  
# fare:female  0.01     0.01  
# ---
#   n = 714, k = 7
bbmle::AICtab(m, m2) # exercise
#    dAIC df
# m2 0    7 
# m  7    4 
sjPlot::plot_model(m2, type = "est")
```
![优势比](R/Rplot122.jpeg)

在上述模型中，大多数系数看起来都很小。为什么？
接着可以进行数据标准化：
现在，我们将用已缩放（除以2个标准差）和居中（减去它们的平均值）的预测值版本重新调整上述模型。对于二进制预测器（`female`），变量将以其平均值为中心，但不按比例缩放。该标准化程序将使系数的大小近似可比。
```r
d$age_scaled <- arm::rescale(d$age)
d$fare_scaled <- arm::rescale(d$fare)
d$female_centered <- arm::rescale(d$female)
# or
# d$female_centered <- d$female - mean(d$female) # same thing

# or:
# m3 <- arm::standardize(m2)
# but we will use arm::rescale so it is clear what we are doing 
```
现在用所有双向交互重新调整模型，但使用年龄和票价的缩放版本：`age_scaled`和`fare_scaled`。这次使用“female”的0-1版本：
```r
m3 <- glm(survived ~ 
    (age_scaled + fare_scaled + female)^2, data = d, family = binomial()) # exercise
arm::display(m3)
# glm(formula = survived ~ (age_scaled + fare_scaled + female)^2, 
#     family = binomial(), data = d)
#                        coef.est coef.se
# (Intercept)            -1.33     0.12  
# age_scaled             -0.72     0.25  
# fare_scaled             0.95     0.32  
# female                  2.51     0.21  
# age_scaled:fare_scaled  0.73     0.60  
# age_scaled:female       1.04     0.41  
# fare_scaled:female      1.40     0.77  
# ---
#   n = 714, k = 7
#   residual deviance = 703.0, null deviance = 964.5 (difference = 261.5)
sjPlot::plot_model(m3, type = "est")
```
![优势比](R/Rplot123.jpeg)

现在拟合同样的模型，但使用`female_centered`：
```r
m4 <- glm(survived ~ 
    (age_scaled + fare_scaled + female_centered)^2, data = d, family = binomial()) # exercise
arm::display(m4)
# glm(formula = survived ~ (age_scaled + fare_scaled + female_centered)^2, 
#     family = binomial(), data = d)
#                             coef.est coef.se
# (Intercept)                 -0.42     0.10  
# age_scaled                  -0.34     0.20  
# fare_scaled                  1.46     0.33  
# female_centered              2.51     0.21  
# age_scaled:fare_scaled       0.73     0.60  
# age_scaled:female_centered   1.04     0.41  
# fare_scaled:female_centered  1.40     0.77  
# ---
#   n = 714, k = 7
#   residual deviance = 703.0, null deviance = 964.5 (difference = 261.5)
sjPlot::plot_model(m4, type = "est")
```
![优势比](R/Rplot124.jpeg)

解读标准化模型
如果泰坦尼克号上的一个人花了250美元买票，那么与另一个花了150美元的人相比，他活下来的几率有多大？（请注意，其中一个标准化预测值的单位变化表示原始变量的2个标准偏差。`fare`的两个标准差约为100美元）。即：
```r
round(sd(d$fare) * 2, 1)
# [1] 105.8
```
我们如何使用“m4”中的系数来确定为一个男人和一个女人的票多支付大约100美元的效果？
请记住，这些模型所做的预测是相同的，只是它们的参数化略有不同。
```r
exp(coef(m4)[["fare_scaled"]] + coef(m4)[["fare_scaled:female_centered"]] * 0.63) # exercise
# [1] 10.38939
exp(coef(m4)[["fare_scaled"]] + coef(m4)[["fare_scaled:female_centered"]] * -0.37) # exercise
# [1] 2.564542
```