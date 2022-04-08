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
			fprintf(fp, "%.4f,", mat[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
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
#endif // DEBUG

	return false;
}
