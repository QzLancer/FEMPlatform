#pragma once

#include <iostream>

/// @brief ���double���͵ľ���
/// @param fname �ļ�����
/// @param rowsize ����
/// @param colsize ����
/// @param mat �����Ա
/// @return �Ƿ�ɹ����
bool printdoubleMatrix(const char* fname, int rowsize, int colsize, double** mat);

/// @brief ���int���͵ľ���
/// @param fname �ļ�����
/// @param rowsize ����
/// @param colsize ����
/// @param mat �����Ա
/// @return �Ƿ�ɹ����
bool printintMatrix(const char* fname, int rowsize, int colsize, int** mat);

/// @brief ���double���͵�����
/// @param fname �ļ�����
/// @param size ��������
/// @param vec ����Ԫ��
/// @return �Ƿ�ɹ����
bool printintVector(const char* fname, int size, int* vec);

/// @brief ���int���͵�����
/// @param fname �ļ�����
/// @param size ��������
/// @param vec ����Ԫ��
/// @return �Ƿ�ɹ����
bool printdoubleVector(const char* fname, int size, double* vec);
