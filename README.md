
# FEMPlatform

#### 介绍
    用于FEM算法研究的平台，纯C++。主要目的是解决针对不同模型求解、不同物理场求解、不同算法研究时，重复造轮子的问题。
    目前主要是实现了NR算法和NDDR算法。

#### 环境需求
	1. visual studio2019 + Windows SDK 10.0以上
	2. Qt 5.14.2 + Qt VS Addin
	3. 适用于最新v142生成工具的C++ MFC(x86和x64)

#### 依赖的第三方库（已经添加到项目中）
    1. Armadillo
    2. SuperLU_MT
    3. Gmsh

#### 构建方法
	打开FEMPlatform_Qt文件夹下的vcxproj文件即可编译。




