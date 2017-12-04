# Matlab常用工具总结

## 二维绘图
### plot函数
#### 函数格式

```
plot(x,y)
```
plot函数以向量x、y为轴，绘制曲线。其中x和y为长度相同的向量。

另有x,y为矩阵的用法，见help。

含多个输入参数的plot函数调用格式为：
```
plot(x1,y1,x2,y2,…,xn,yn)
```
- 当输入参数都为向量时，x1和y1，x2和y2，…，xn和yn分别组成一组向量对，每一组向量对的长度可以不同。每一向量对可以绘制出一条曲线，这样可以在同一坐标内绘制出多条曲线。
- 当输入参数有矩阵形式时，配对的x,y按对应列元素为横、纵坐标分别绘制曲线，曲线条数等于矩阵的列数。

#### 图形保持 hold
在一个figure上面先后plot多次。
1. hold on:启动图形保持功能，当前坐标轴和图形都将保持，此后绘制的图形都将添加在这个图形之上，并且自动调整坐标轴的范围。
2. hold off:关闭图形保持功能。
3. hold ：在hold on 和hold off命令之间进行切换。

#### 设置曲线样式格式
详见**help LineSpec**

用于确定所绘曲线的线型、颜色和数据点标记符号。

当选项缺省时，MATLAB规定，线型一律用实线，颜色将根据曲线的先后顺序依次使用。


```
plot(x1,y1,’cms’,....)
```
c为颜色，m为标记，s为线型。
```
plot(x,y,'-.om')
%'-.'line,'o'circle marker,'m'magenta color 洋红色
```
另可设置Linewidth,MarkerEdgeColor,MarkerFaceColor,MarkerSize

```
figure
plot(t,sin(2*t),'-mo','Linewidth',2,'MarkerEdgeColor', 'k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10)
```
当然，亦可于之后通过绘图对话窗的GUI界面对上述参数进行修改。



#### 图形标记
在绘制图形的同时，可以对图形加上一些说明，如图形名称、图形某一部分的含义、坐标说明等。

```
title(‘加图形标题’);
%当前轴的正上方居中位置处输出文本作为标题    
xlabel('加X轴标记');     
ylabel('加Y轴标记');       
text(X,Y,'添加文本');
```
#### 坐标axis、网格grid、边框box控制
axis函数的调用格式为：

```
axis([xmin xmax ymin ymax zmin zmax])
```
axis函数功能丰富，常用的格式还有：
纵、横坐标轴采用等长刻度:
```
axis equal
```
产生正方形坐标系(缺省为矩形)。
```
axis square
```
使用缺省设置。
```
axis auto
```
取消坐标轴。
```
axis off
```
显示坐标轴。
```
axis on
```
相似的还有
- 网格控制grid：grid on/off命令控制是画还是不画网格线，不带参数的grid命令在两种状态之间进行切换。
- 边框控制box：box on/off命令控制是加还是不加边框线，不带参数的box命令在两种状态之间进行切换。

#### 图例 legend

```
x=0:pi/100:2*pi;
y1=sin(x);
y2=cos(x);
plot(x,y1,x,y2, '--');
legend('sin(x)','cos(x)');
```
#### 极坐标图 polar
调用格式
```
polar(theta,rho,选项)
```
## 三维绘图
### 三维网格图 mesh

```
mesh(x,y,z,'LineSpec')
```
示例

```
x=linspace(-2, 2, 25); % 在x轴上取25点
y=linspace(-2, 2, 25); % 在y轴上取25点
[xx,yy]=meshgrid(x,y); % 形成网格，xx和yy都是25x25的矩阵
zz=xx.*exp(-xx.^2-yy.^2);
%计算函数值，zz也是25x25的矩阵
mesh(xx, yy, zz); %画出立体网状图
```
### 三维曲面图 surf
各线条之间的补面用颜色填充。


```
surf(x,y,z,'LineSpec')
```
## 灰度图 imagesc,colormap(gray)
三维向二维投影
```
figure
imagesc(xpart,ypart,xy_data(:,:,z));
colormap(gray(64))
colorbar
```
关于colorbar 的一些常用设置：
1. 在colorbar所在figure中，选择文件->导出设置->字体，即可根据要求修改字体大小。
2. 对于colorbar的范围设置，需要使用语句
caxis([m,n])，m,n即为需要的范围。
3. 当然colorbar可以隐去数值，只保留颜色，使用语句
colorbar('YTickLabel',{' '})
即可。
4.翻转xy轴使用 'axis xy'



## 动画 getframe,moviein,movie

```
%动画功能函数：getframe、moviein和movie
%播放一个不断变化的眼球程序。
m=moviein(20); %建立一个20个列向量组成的矩阵
for j=1:20
  plot(fft(eye(j+10)))
  %绘制出每一幅眼球图并保存到m矩阵中
  m(:,j)=getframe;
end
movie(m,10);%以每秒10幅的速度播放画面
```


## 其他
### cell模式
在一个长长的脚本m文件中，可能需要对其中的一段反复修改，查看执行效果，这时，cell模式就非常有用了。

cell模式相当于将其中的代码拷贝到命令窗口中运行。两个%后接一个空格` %% `开始一个cell。

将输入光标放到一个cell中时，背景将变为浅黄色，` Ctrl+Enter`执行cell中的代码。

执行cell中代码时不需要保存m文件，该m文件可以不在路径列表中。

cell模式中，断点不起作用，当然，调用的子程序中的断点还是正常的。

### 随机数生成函数
` rand(m,n)`生成均匀分布的 m * n 矩阵
注： ` rand(n,n)=rand(n,n)`

` randn(m,n)` 满足正态分布的m * n 矩阵

` randperm(m)` 生成由1:m组成的随机序列

` perms(1:n)` 生成由1:n组成的全排列，共n！个。

` random(`name`,A1,A2,A3,m,n)` 满足name分布的m * n随机矩阵，A1...A3等是所需的参数
name可以取值为
- norm or Normal
- unif or Uniform
- beta or Beta
- exp or Exponential
- gam or Gamma
- geo or Geometric
- unid or Discrete Uniform

具体参见help random

### 取整函数
` fix(x)` 截尾取整，舍掉小数部分
` floor(x)` 不超过x的最大整数
` ceil(x)` 不小于x的最小整数
` round(x)` 四舍五入取整

### 转置
对矩阵A，`A.'`为普通转置，` A'`为共轭转置。

当然对于是实数阵无差别。

### 创建行向量(row vector)
#### 已知间隔
指定步长s
```
x = 0:pi/100:2*pi
y = [0:0.5:360]*pi/180
```
#### 已知所需元素数
linspace（to generate linearly spaced vectors）

```
y = linspace(a,b,n)
% generate n points linearly spaced between and including a and b
```
### save/load
#### save
有命令模式和函数模式（即加括号和单引号）两种
- 命令模式
```
save tian.mat x y
```
- 函数模式
```
save('tian.mat','x','y')
```
#### load
用load函数，可以将数据读入到matlab的工作空间中。

```
load('tian.mat')
```

注意，load 是将filename.mat中的所有变量读入matlab工作空间中,也可以选择读入哪个变量。注意原mat保存的哪一个（些）变量，**变量名字是没有变的**。

还有一个问题是我的mat文件中保存这一个变量,可是默认的读入matlab中后,还是保存时用的名字,但是,我想用一个新名字代替,怎么办?具体用程序描述如下:

```
 save tian.mat p
```

那么load tian.mat之后，就可以在工作空间中看到p变量了。可是，如果用load读入之后，不想用变量名p了，怎么办？

具体解决办法：

```
s=load('tian.mat');

sc=struct2cell(s);

t=cell2mat(sc);
```


那么,读入的struct类型变量就被转换成cell类型数据,然后再转换为double类型的数据。

---

参考文档
1. Matlab help browser.
2. [Matlab绘图命令](http://www.cnblogs.com/hxsyl/archive/2012/10/10/2718380.html)
3. [Matlab中save与load函数的使用](http://www.cnblogs.com/rong86/p/3559861.html)
