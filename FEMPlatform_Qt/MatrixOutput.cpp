#include "MatrixOutput.h"

bool printdoubleMatrix(const char* fname, int rowsize, int colsize, double** mat)
{
#ifdef DEBUG
	FILE* fp;
	char fullpath[256];
	sprintf(fullpath, "../matrix/%s", fname);
	fp = fopen(fullpath, "w+");
	for (int i = 0; i < rowsize; ++i)
	{
		for (int j = 0; j < colsize; ++j) {
			fprintf(fp, "%.12f,", mat[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return true;
#endif // DEBUG

	return false;
}

bool printintMatrix(const char* fname, int rowsize, int colsize, int** mat)
{
#ifdef DEBUG
	FILE* fp;
	char fullpath[256];
	sprintf(fullpath, "../matrix/%s", fname);
	fp = fopen(fullpath, "w+");
	for (int i = 0; i < rowsize; ++i)
	{
		for (int j = 0; j < colsize; ++j) {
			fprintf(fp, "%d,", mat[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return true;
#endif // DEBUG

	return false;
}

bool printintVector(const char* fname, int size, int* vec)
{
#ifdef DEBUG
	FILE* fp;
	char fullpath[256];
	sprintf(fullpath, "../matrix/%s", fname);
	fp = fopen(fullpath, "w+");
	for (int i = 0; i < size; ++i)
	{
		fprintf(fp, "%d\n", vec[i]);
	}
	fclose(fp);
	return true;
#endif // DEBUG
	return false;
}

bool printdoubleVector(const char* fname, int size, double* vec)
{
#ifdef DEBUG
	FILE* fp;
	char fullpath[256];
	sprintf(fullpath, "../matrix/%s", fname);
	fp = fopen(fullpath, "w+");
	for (int i = 0; i < size; ++i)
	{
		fprintf(fp, "%.12f\n", vec[i]);
	}
	fclose(fp);
	return true;
#endif // DEBUG

	return false;
}
