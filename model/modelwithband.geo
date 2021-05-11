lc = 1e-3;
d = 0.0004;
//衔铁位置，从0到6，位置为从下到上，气隙0.0025，最终位置行程0.0024，相距0.0001mm
n = 0;

//铁磁部分几何点
Point(1) = {0.00265, 0.0099, 0, lc};
Point(2) = {0.0235, 0.0099, 0, lc/5};
Point(3) = {0.0235, 0.011, 0, lc/5};
Point(4) = {0.02395, 0.011, 0, lc/5};
Point(5) = {0.00265, 0.0077, 0, lc/11};
Point(6) = {0.004735, 0.0077, 0, lc/20};
Point(7) = {0.004735, 0.0083, 0, lc/5};
Point(8) = {0.0235, 0.0083, 0, lc/5};
Point(9) = {0.0095565, 0.00685, 0, lc};
Point(10) = {0.0212435, 0.00685, 0, lc};
Point(11) = {0.0095565, -0.00785, 0, lc};
Point(12) = {0.0212435, -0.00785, 0, lc};
Point(13) = {0.0024, 0.0052+n*d, 0, lc/11};
Point(14) = {0.004735, 0.0052+n*d, 0, lc/20};
Point(15) = {0.004835, 0.0027, 0, lc/20};
Point(16) = {0.00595, 0.0027, 0, lc/2};
Point(17) = {0.00665, 0.0015, 0, lc/2};
Point(18) = {0.0016, 0.0001+n*d, 0, lc/2};
Point(19) = {0.0024, 0.0001+n*d, 0, lc};
Point(20) = {0.0016, -0.0088+n*d, 0, lc/10};
Point(21) = {0.004735, -0.0088+n*d, 0, lc/20};
Point(22) = {0.00665, -0.009, 0, lc};
Point(23) = {0.019, -0.009, 0, lc/4};
Point(24) = {0.019, -0.01, 0, lc/4};
Point(25) = {0.0235, -0.01, 0, lc/4};
Point(26) = {0.004835, -0.0108, 0, lc/10};
Point(27) = {0.02395, -0.0108, 0, lc/4};
Point(28) = {0.00595, 0.0015, 0, lc/5};

//圆心
Point(29) = {0, 0, 0, lc};

//内圈边界
Point(30) = {0, 0.04, 0, 5*lc};
Point(31) = {0, -0.04, 0, 5*lc};
Point(32) = {0.04, 0, 0, 5*lc};

//外圈边界
Point(33) = {0, 0.05, 0, 5*lc};
Point(34) = {0, -0.05, 0, 5*lc};
Point(35) = {0.05, 0, 0, 5*lc};

//重分网空气区域
Point(36) = {0 ,0.0077 , 0, lc/2};
Point(37) = {0.004835,0.0077 , 0, lc/20};
Point(38) = {0, -0.009, 0, (n+1)*lc/10};
Point(39) = {0.004835, -0.009,0, lc/20};


//构建线条
Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,27} ;
Line(5) = {1,5} ;
Line(6) = {5,6} ;
Line(7) = {6,7} ;
Line(8) = {7,8} ;
Line(9) = {8,25} ;
Line(10) = {9,10} ;
Line(11) = {9,11} ;
Line(12) = {11,12} ;
Line(13) = {10,12} ;

//衔铁区域线条
Line(14) = {13,14} ;
Line(15) = {13,19} ;
Line(16) = {18,19} ;
Line(17) = {18,20} ;
Line(18) = {20,21} ;
Line(19) = {14,21} ;

//轭铁区域线条
Line(20) = {15,16} ;
Line(21) = {16,28} ;
Line(22) = {28,17} ;
Line(23) = {39,26} ;
Line(24) = {17,22} ;
Line(25) = {22,23} ;
Line(26) = {23,24} ;
Line(27) = {24,25} ;
Line(28) = {26,27} ;

//内圈圆边1
Circle (29) = {30, 29, 32} ;
Circle (30) = {32, 29, 31} ;


//外圈圆边2
Circle (31) = {33, 29, 35} ;
Circle (32) = {35, 29, 34} ;

//轴对称边界
Line(33) = {33,30};
Line(34) = {30,36};
Line(35) = {38,31};
Line(36) = {31,34};

//重分网空气区域边界
Line(37) = {36, 5};
Line(38) = {6, 37};
Line(39) = {37, 15};
Line(40) = {38, 39};
Line(41) = {36, 29};
Line(42) = {29, 38};
Line(43) = {15, 39};

//Curve Loop Define
//衔铁区域
Curve Loop(1) = {14, 19, -18, -17, 16, -15} ;
Plane Surface(1) = {1} ;

//轭铁区域
Curve Loop(2) = {-5, 1, 2, 3, 4, -28, -23, -43, 20, 21, 22, 24, 25, 26, 27, -9, -8, -7, -6} ;
Plane Surface(2) = {2} ;

//线圈区域
Curve Loop(3) = {10, 13, -12, -11} ;
Plane Surface(3) = {3} ;

//内部空气域
Curve Loop(4) = {-10, -13, 12, 11} ;
Curve Loop(5) = {-38, 7, 8, 9, -27, -26, -25, -24, -22, -21, -20, -39};
Plane Surface(4) = {4,5} ;

//重分网空气域
Curve Loop(6) = {14, 19, -18, -17, 16, -15} ;
Curve Loop(7) = {-41, 37, 6, 38, 39, 43, -40, -42};
Plane Surface(5) = {6,7} ;

//空气圆环
Curve Loop(8) = {-34, 29, 30, -35, 40, 23, 28, -4, -3, -2, -1, 5, -37} ;
Plane Surface(6) = {8};

//半圆内部空气区域
Curve Loop(9) = {-33, 31, 32, -36, -30, -29};
Plane Surface(7) = {9};





