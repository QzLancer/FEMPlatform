#pragma once

#include <iostream>

/// @brief 输出double类型的矩阵
/// @param fname 文件名称
/// @param rowsize 行数
/// @param colsize 列数
/// @param mat 矩阵成员
/// @return 是否成功输出
bool printdoubleMatrix(const char* fname, int rowsize, int colsize, double** mat);

/// @brief 输出int类型的矩阵
/// @param fname 文件名称
/// @param rowsize 行数
/// @param colsize 列数
/// @param mat 矩阵成员
/// @return 是否成功输出
bool printintMatrix(const char* fname, int rowsize, int colsize, int** mat);

/// @brief 输出double类型的向量
/// @param fname 文件名称
/// @param size 向量长度
/// @param vec 向量元素
/// @return 是否成功输出
bool printintVector(const char* fname, int size, int* vec);

/// @brief 输出int类型的向量
/// @param fname 文件名称
/// @param size 向量长度
/// @param vec 向量元素
/// @return 是否成功输出
bool printdoubleVector(const char* fname, int size, double* vec);
