#include "SluMTMatrixSolver.h"


SluMTMatrixSolver::SluMTMatrixSolver() :
    nprocs(std::thread::hardware_concurrency())
{
}

SluMTMatrixSolver::~SluMTMatrixSolver()
{
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    Destroy_SuperMatrix_Store(&A);
}

double* SluMTMatrixSolver::solveMatrix(vector<vector<int>> locs, vector<double> vals, vector<double> F, int valsize, int vecsize)
{
    printf("nprocs: %d\n", nprocs);
    if (L.ncol > 0 || L.nrow > 0) {
        Destroy_SuperNode_Matrix(&L);
    }
    if (U.ncol > 0 || U.nrow > 0) {
        Destroy_CompCol_Matrix(&U);
    }

    int m, n, nnz, info, panel_size = 0;
    int* xa = nullptr, * asub = nullptr, * perm_r = nullptr, * perm_c = nullptr;
    double* a = nullptr, * rhs = nullptr, * res = nullptr;
    superlu_memusage_t superlu_memusage;
    //将数据转化成列压缩存储形式
    m = vecsize; n = vecsize;
    std::map<int, std::map<int, double>> mapmatrix; //mapmatrix[列编号][行编号][值]
    int row, col;
    double val;
    for (int i = 0; i < valsize; ++i) {
        row = locs[0][i];
        col = locs[1][i];
        val = vals[i];
        if (mapmatrix.count(col) == 0) {
            std::map<int, double> temp;
            temp[row] = val;
            mapmatrix[col] = temp;
        }
        else {
            if (mapmatrix[col].count(row) == 0) {
                mapmatrix[col][row] = val;
            }
            else {
                mapmatrix[col][row] += val;
            }
        }
    }

    nnz = 0;
    xa = intMalloc(vecsize + 1);
    xa[0] = 0;
    for (std::map<int, std::map<int, double>>::iterator m = mapmatrix.begin(); m != mapmatrix.end(); ++m) {
        nnz += m->second.size();
        xa[(m->first) + 1] = nnz;
    }
    asub = intMalloc(nnz);
    a = doubleMalloc(nnz);
    int i = 0;
    for (std::map<int, std::map<int, double>>::iterator m = mapmatrix.begin(); m != mapmatrix.end(); ++m) {
        for (std::map<int, double>::iterator iter = m->second.begin(); iter != m->second.end(); ++iter) {
            asub[i] = iter->first;
            a[i] = iter->second;
            //            cout << "a:" << i << " = " << a[i] << endl;
            ++i;
        }
    }

    dCreate_CompCol_Matrix(&A, n, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    int nrhs = 1;
    if (!(rhs = doubleMalloc(m * nrhs))) printf("Malloc fails for rhs[].\n");
    for (i = 0; i < m; ++i) rhs[i] = F[i];
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_r = intMalloc(m))) printf("Malloc fails for perm_r[].\n");
    if (!(perm_c = intMalloc(n))) printf("Malloc fails for perm_c[].\n");

    int permc_spec = 1;
    get_perm_c(permc_spec, &A, perm_c);

    //n = A.ncol;
    pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);

    //读取结果
    if (info == 0) {
        res = doubleMalloc(m);
        for (int i = 0; i < vecsize; ++i) {
            res[i] = ((double*)((DNformat*)B.Store)->nzval)[i];
            //printf("res %d : %f\n", i, res[i]);
        }
    }
    else {
        printf("dgssv() error returns INFO= " IFMT "\n", info);
        if (info <= n) { /* factorization completes */
            superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
            printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
                superlu_memusage.for_lu / 1e6,
                superlu_memusage.total_needed / 1e6,
                superlu_memusage.expansions);
        }
        exit(0);
    }

    SUPERLU_FREE(perm_c);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(rhs);
    SUPERLU_FREE(a);
    SUPERLU_FREE(asub);
    SUPERLU_FREE(xa);

    return res;
}

void SluMTMatrixSolver::setNumberofProcs(const int _nprocs)
{
    nprocs = _nprocs;
}

