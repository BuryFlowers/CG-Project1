# CG Project1


## 项目简介

本项目为计算机图形学秋学期大作业，算法方面分别实现了扫描线zbuffer， 层次zbuffer（简单模式）以及层次zbuffer（完整模式）。使用到了OpenGL的glm库，编写预览的可编程管线时参考了[LearnOpenGL](https://learnopengl-cn.github.io/)的代码，保存离线渲染结果图片时使用了[tooJpeg](https://create.stephan-brumme.com/toojpeg/)的库，其余所有代码均在根据图形学课件的算法描述由本人手写。


## 程序介绍

由于上述三种算法全都用CPU实现，在渲染大型模型时无法达到实时渲染的要求，所以本人决定使用OpenGL的可渲染管线对模型进行预览，在选择好合适的相机位置与角度后，再使用三种算法进行离线渲染。
在程序加载模型完毕后，预览窗口会自动出现，使用wasd可以控制相机位置，使用方向键上下可以控制相机移动速度，使用鼠标滑动可以控制相机角度，使用鼠标滚轮可以放大或者缩小画面。
按下“P”键可以输出相机具体参数，按下“R”键可以开始渲染当前画面，在渲染结束后会输出三种算法各自的用时和剔除的三角形数量，在result文件夹可以看到三种算法的渲染结果。


## 编写环境

本程序使用Microsoft Visual Studio Community 2022 (64 位) 编写。
CPU为Intel(R) Core(TM) i7-7700HQ。
内存16.0GB。

## 算法描述

为了简化问题，本项目只支持渲染三角形面片。

### 扫描线zbuffer
本算法在窗口空间进行扫描，对于每条扫描线，记录当前与扫描线相交的所有三角形的两边，从左端点到右端点更新zbuffer与颜色缓冲。

#### 时间复杂度
本质上与普通的zbuffer算法没有数量级上的差别

#### 优点
在渲染的面片的互相遮挡程度不高时效率很高。

#### 缺点
数据结构与精度问题复杂。每条与扫描线相交的三角形的两边对当前扫描线的zbuffer进行了多次更新，仍有优化空间。

### 层次zbuffer（简单模式）
在普通zbuffer的基础上建立层次zbuffer，高一层zbufer的值为相邻四个低层非空zbuffer的最大值（最远值）。当检查一个物体的可见度时，从最高层zbuffer进行检查，若当前层的zbuffer内的所有像素都被渲染，且当前zbuffer的值小于物体的最小z值，则代表该物体必定被遮挡，可以提前剔除。

#### 时间复杂度
若被遮挡的三角形较多，则渲染效率较高。

#### 优点
运用三角形的遮挡关系，可以在渲染前提前剔除某些三角形。

#### 缺点
效率与渲染三角形的顺序较为相关，若按z值顺序从后往前渲染，则基本剔除不了多少三角形。且依赖于物体在屏幕中的占比，若物体面片非常多，且全部聚集在很小的一部分像素中，那么层次zbuffer会退化成普通的zbuffer。

### 层次zbuffer（完整模式）
在简单模式的基础上建立八叉树，按照空间相对位置将三角形进行分类，提前渲染离相机较近的与上一帧渲染没有被剔除的三角形，若当前八叉树节点被层次zbuffer拒绝，则该节点下所有三角形都可以被直接剔除。

#### 时间复杂度
效率与被剔除的三角形数目挂钩。

#### 优点
与简单模式相比，容易剔除更多的三角形。

#### 缺点
并没有解决层次zbuffer在物体占比较小时的退化问题。


## 效率测试
本次测试分别渲染了拥有200/5k/100k/2m个三角形的模型，分别从远近2个位置进行渲染，每次渲染5张720P的图片，记录平均渲染时间以及加速比，与剔除的面片数量与其占比。

### 200面片
第一个测试使用一个具有200个三角形面片的兔子模型，其近距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/200/near/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/200/near/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/200/near/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 4.56ms | 60.71ms (+1231%) | 61.71ms (+1253%) |
| 剔除三角面数量 | 0 | 140 (70%) | 140 (70%) |

远距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/200/far/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/200/far/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/200/far/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 2.40ms | 43.23ms (+1701%) | 34.71ms (+1346%) |
| 剔除三角面数量 | 0 | 0 | 4 (2%) |

### 5k面片
第二个测试使用一个具有5k个三角形面片的兔子模型，其近距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/5k/near/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/5k/near/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/5k/near/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 12.19ms | 74.55ms (+511%) | 81.90ms (+572%) |
| 剔除三角面数量 | 0 | 3900 (78%) | 3915 (78%) |

远距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/5k/far/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/5k/far/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/5k/far/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 5.88ms | 21.71ms (+269%) | 24.02ms (+309%) |
| 剔除三角面数量 | 0 | 2410 (48%) | 2285 (46%) |

### 100k面片
第三个测试使用一个具有100k个三角形面片的相机模型，其近距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/100k/near/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/100k/near/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/100k/near/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 88.14ms | 139.36ms (+58%) | 133.21ms (+51%) |
| 剔除三角面数量 | 0 | 47337 (47%) | 72118 (72%) |

远距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/100k/far/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/100k/far/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/100k/far/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 71.83ms | 77.12ms (+7%) | 89.03ms (+24%) |
| 剔除三角面数量 | 0 | 26859 (27%) | 54750 (55%) |

### 2m面片
第四个测试使用一个具有2m个三角形面片的雕塑模型，其近距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/2m/near/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/2m/near/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/2m/near/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 1,323.80ms | 704.24ms (-47%) | 606.22ms (-54%) |
| 剔除三角面数量 | 0 | 814763 (41%) | 1392303 (70%) |

远距离的渲染结果如下表所示：

| | 扫描线zbuffer | 层次zbuffer（简单模式） | 层次zbuffer（完整模式） |
| ----------- | ----------- | ----------- | ----------- |
| 渲染图片 | ![扫描线zbuffer](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/2m/far/scanline_zbuffer.jpg) | ![层次zbuffer（简单模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/2m/far/simple_hierachical_zbuffer.jpg) | ![层次zbuffer（完整模式）](https://github.com/BuryFlowers/CG-Project1/blob/master/CG%20Project1/images/2m/far/hierachical_zbuffer.jpg) |
| 渲染平均用时 | 1,012.91ms | 564.56ms (-44%) | 462.46ms (-54%) |
| 剔除三角面数量 | 0 | 1050100 (53%) | 1580420 (79%) |

## 结果分析
当面片数量较少时，层次zbuffer由于需要额外的处理，所以效率远不如扫描线zbuffer。而当面片数量很多时，层次zbuffer就会快于扫描线zbuffer，且完整模式一般情况下也会快于简单模式，与扫描线zbuffer相比，平均用时最快可以减少50%~60%。然而，当有许多面片聚集在少数像素中时，层次zbuffer的加速效果便没有那么明显了。
