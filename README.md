
# FEMPlatform

## 介绍
用于FEM算法研究的平台，纯C++。主要目的是解决针对不同模型求解、不同物理场求解、不同算法研究时，重复造轮子的问题。
目前主要是实现了NR算法、NDDR算法和Schwarz算法，包含二维平面和二维轴对称低频电磁的计算。

## 环境需求
1. visual studio2019及以上：https://visualstudio.microsoft.com/zh-hans/
2. Windows SDK 10.0以上+MSVC v142生成工具+适用于最新v142生成工具的C++ MFC(x86和x64)：打开visual studio installer, 安装如下组件
![](README_md_files/5816ae10-2384-11ed-98d1-41352cacc72b.jpeg?v=1&type=image)
3. Qt 5.14.2 + Qt VS Addin：https://download.qt.io/archive/qt/
Qt安装时选择MSVC2019 64bit编译器
5. CUDA v11.1及以上版本（支持shared memory）+Visual studio拓展插件

## 依赖的第三方库（已经添加到项目中，在win_x64下不需要重新编译）
1. Armadillo（矩阵运算库）
2. SuperLU_MT（并行直接法矩阵求解库）
3. Gmsh（网格剖分库）
4. Metis（区域分解库）

## 构建方法
打开FEMPlatform_Qt文件夹下的vcxproj文件即可编译。（FEM_Platform为早期版本，FEMPlatform_Qt为使用Qt进行构建的后期版本，导入了Gmsh中的重分网模块，以实现基于重分网技术的电磁动态特性计算）
	
## 代码的主要类

### FEMCore

模块与`main()`函数之间的接口，设计`FEMCore`类的目的主要是为了后续图形化界面的开发，在`main()`函数中通过调用`FEMCore`的成员配置求解参数，并调用`FEMCore::solve()`求解FEM问题。

![](README_md_files/da620480-2386-11ed-8593-e54b1b5abc89.jpeg?v=1&type=image)



### FEMModel
该类实现模型的构建，需要创建新模型时，创建一个子类继承`FEMModel`类，然后设置模型的全部必要参数。此处采用template设计模式，基类提供虚函数用于设置参数，并且在`init()`函数中统一调用这些虚函数进行初始化。这样，一方面基类作为抽象类，无法直接构建基类的实例；另一方面，继承自该类的其它类就必须override父类的虚函数，并最终调用`init()`函数，避免了未初始化参数的存在。`FEMModel`这一基类还提供一些获取几何、网格、激励等数据的接口。

接下来通过案例简述模型的初始化过程

`FEMContactorNonLinearModel`类构建了二维轴对称电磁静态特性计算模型，该类继承自`FEMModel`类，并且对FEMModel中的纯虚函数进行重载。各个函数的具体实现参考`FEMContactorNonLinearModel.c`。由于静态特性计算不需要设置运动区域和相关的弹簧条件，因此`buildGeometry2Deformed()`和`buildGeometry2MovingPart()`的设置不影响计算，可以直接return。
![](README_md_files/4a89c000-2354-11ed-98d1-41352cacc72b.jpeg?v=1&type=image)

`RelayDynamicModel`类构建了二维轴对称电磁机构瞬态特性计算模型，相比于静态特性计算，瞬态特性的计算需要额外对电路参数和弹簧特性进行初始化，并指定空气域和运动区域。电路相关参数在`createElement2Material()`中进行初始化。

模型构建过程还涉及到如下关键类：
 - `FEMBoundary`：边界条件类，目前包含轴对称边界条件和自然边界条件。
 - `FEMCoil`：线圈结构体，用于储存线圈关键参数。
 - `FEMMovingPart`：运动部件参数存储及弹簧力计算。
 - `FEMMaterial`：保存模型中材料的基本信息，包含材料名称、材料的B-H曲线、线圈参数等。除此之外，材料管理模块还负责材料基本参数计算，并提供接口，例如计算磁导率对磁感应强度的偏导。

类中成员函数及变量在此不赘述，在该类的相关代码中给出了详细备注。

### FEMMeshManager

读取COMSOL生成的`.mphtxt`网格文件，或读取Gmsh生成的`.msh`网格文件，或调用Gmsh接口对`.geo`几何文件分网，并将网格数据文件解析成可用于有限元分析的数据结构。针对电磁机构动态特性的计算，`FEMMeshManager`还负责保存模型的重分网单元信息以及衔铁单元信息，并通过这些信息在模块中实现模型的局部重分网。

网格剖分还涉及到如下关键类：
- `CNode`：用于节点信息的存储，整合了NDDR算法所需的信息
- `CVtxElement, CEdgElement, CTriElement`：节点单元、线单元和三角形单元信息

上述类存储在`FEMDataType.h`中
- `Gmsh`：负责.geo几何文件分网和重分网，其核心函数通过静态成员的形式进行封装。只有Qt版本的代码导入了Gmsh模块。

派生`FEM2DMeshManager`和`FEM3DMeshManager`，目前只用到`FEM2DMeshManager`，实现二维网格的剖分和网格文件读取。

### FEMSolver

代码的核心部分，用于实现有限元求解流程，派生`FEM2DSolver`和`FEM3DSolver`，再进一步派生`FEM2DNRSolver, FEM2DNDDRSolver, FEM2DNDDRCUDASolver, FEM2DSchwarzSolver`。每个类中包含静态特性计算和动态特性计算的函数

`FEM2DNRSolver`：经典的全局NR迭代，没什么难度

`FEM2DNDDRSolver`：CPU版本的NDDR求解器，`FEM2DNDDRSolver::solveStatic()`函数中，针对NDDR的静态特性计算，提供了很多实现的版本，但是收敛性未能得到改善，供参考。

`FEM2DNDDRCUDASolver`：GPU版本的NDDR求解器，使用GPU求解NDDR时需要注意，要将内存中的变量完全拷贝至GPU中，通过_global_关键字指定使用GPU并行执行的函数，并且调用的子函数必须是_device_版本。同样收敛性差。

`FEM2DSchwarz`：Schwarz重叠子域区域分解的实现，通过metis实现非重叠区域分解后，进一步得到了一个重叠单元的子域区域分解，尽管实现了并行化并且收敛性好，但是求解效率相比于传统的NR方法几乎没有提升。

求解过程中涉及到的关键类：
- `MatrixSolver`：矩阵求解器，用于[S][A]=[F]线性系统的计算，`MatrixSolver`为抽象类，提供`solveMatrix()`这一接口函数，S输入的是三元组表示法的稀疏矩阵。目前派生`SluMTMatrixSolver`类，通过调用SuperLU_MT这一并行化的直接法矩阵求解器解算方程组。

### MatrixOutput

用于矩阵调试，内置4个函数：`printdoubleMatrix(),  printintMatrix(), printintVector(), printdoubleVector()`。输出的是.csv格式文件，可以直接通过excel或其它csv查看工具查看输出的矩阵。

## 代码的主要功能实现

### 用于有限元计算的网格数据的解析

有限元分析需要如下网格数据：

 - 节点编号及对应的几何坐标
 - 单元中包含的的节点编号，当前代码版本中存在点单元（一个节点）、线单元（两个节点）、一阶三角形面单元（三个节点）
 - 单元所在的域的编号，用于加载负载、材料和定位边界条件，例如对于二维模型仿真，线单元所在域的编号可用于定位边界条件，面单元所在域的编号可用于加载负载和材料。

以COMSOL导出的网格文件为例，在COMSOL建模时，可以看到对于几何的边界和连通域，都有一个单独的域编号，网格导出时，各个单元对应的域编号信息也被包含在导出的网格内。因此可以通过解析并保存相关的域单元数据，给不同的域指定不同单元的边界条件、材料及负载等信息。
![](README_md_files/0d3a1210-23ac-11ed-8593-e54b1b5abc89.jpeg?v=1&type=image)
![](README_md_files/1daae700-23ac-11ed-8593-e54b1b5abc89.jpeg?v=1&type=image)

### 核心算法（NR、NDDR，静态特性及瞬态特性计算）
参考硕士毕业论文：https://kns.cnki.net/kcms2/article/abstract?v=3uoqIhG8C475KOm_zrgu4lQARvep2SAkueNJRSNVX-zc5TVHKmDNkn_cgVevTER5o8K00EeWCJqvqWb7Qp25O6Dvqiy495u9&uniplatform=NZKPT

### 基于.vtk的后处理

软件求解完成后输出.vtk文件，可以通过paraview绘制云图、等高线图、矢量图和衔铁的运动动画等。

![](README_md_files/d118e490-23ac-11ed-8593-e54b1b5abc89.jpeg?v=1&type=image)

## 代码存在的问题

 1. 当前在 `FEMCore::solve(string analysistype)`中通过analysistype字符串来判断是计算静态特性还是动态特性，最合理的方法是放在`FEMCore::setAnalysisType(string analysistype)`中判断。
 2. `solve`函数执行过程中，将存储网格信息的数组一个个从`meshmanager`传递至`solver`中，代码冗余，并且为了节约空间和时间，采用的是浅拷贝的方式，存在内存泄漏或指针被多次析构的隐患。可以考虑直接将`meshmanager`传入`solver`中，在`solver`中处理网格信息。
 3. 后处理，当前并不是完全通过`FEMCore`的`postprocess()`函数实现，而是将很多后处理内容放在了`FEMSolver`类中，因为考虑到动态特性计算时，求解迭代步数需要写入文件名称等问题。建议将后处理相关内容单独整理成`PostProcess`类。
 4. 动态特性计算，算法方面，电路-磁场-机械运动耦合过程，当前采用的策略是在每一个时间步计算电路-磁场运动过程，使其收敛，再用后向欧拉法计算机械运动（采用的也不是严格的后向欧拉法），从ode求解的角度应该还有优化的空间。
 5. 动态特性计算时，指定重分网区域、电流加载区域和运动区域，原本的设想是在`FEMModel`指定后，求解过程中直接读取即可，不需要重新指定，但是在动态特性求解过程时，需要在`solveWeakCouple()`函数中重新指定电流所在区域，`solveSpringForce()`函数中重新指定衔铁所在区域，违背了修改模型后用于求解的核心代码不再需要修改的原则，因此相关代码需要进一步优化。
 6. 采用Release模式计算NR算法时，程序运行时报错，猜测是SuperLU_MT Release版本动态库编译的问题。
 7. FEMSolver存在类过多的问题，如果拓展三维模型计算，添加物理场，新增的类将更多，建议将算法抽象到`FEMSolveStrategy`中。
 8. NDDR算法在BH曲线非线性程度高时不收敛，推测原因是Jacobi矩阵不完整，具体参考文章《NDDR方法的收敛性分析》。
 9. GPU版本的NDDR代码中，用于收敛性分析的规约求和对求解效率的影响很大。
