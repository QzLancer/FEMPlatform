#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <Algorithm>
#include <string>

const int CudaThrdNum = 64;
const int CudaBlckNum = 32;
const double CurrentDensity = 8e6;
const int maxitersteps = 20000;
#define Mu0 1.256637e-6;
#define PI 3.1415927;
const double maxerror = 1e-9;

using namespace std;

struct CNode
{
    double x{ 0 }, y{ 0 }, z{ 0 };
    int bdr{ 0 }; //边界条件
    double A{ 0 };
    double A_old{ 0 };
    double At{ 0 };
    double At_old{ 0 };
    int NumberofNeighbourElement{ 0 };
    int NeighbourElementId[15];	//和节点相关的单元编号
    int NeighbourElementNumber[15];	//节点在对应单元中的编号
    double NodeForcex{ 0 }, NodeForcey{ 0 }, NodeForcez{ 0 }, NodeForce{ 0 };
};
struct CVtxElement {
    int n{ 0 };
    int domain{ 0 };
};

struct CEdgElement {
    int n[2]{ 0 };
    double x[2]{ 0 };
    double y[2]{ 0 };
    double z[2]{ 0 };
    int domain{ 0 };
};

struct CTriElement {
    int n[3]{ 0 };// ni, nj, nk;//
    double Q[3]{ 0 };// Qi, Qj, Qk;
    double R[3]{ 0 };// Ri, Rj, Rk;
    double C[3][3];// 单元系数矩阵
    double area{ 0 }; 
    double rc, zc;
    double xdot;
    int domain{ 0 };
    double J{ 0 };  //负载，暂时只考虑电流，直流，后续需要单独建一个类
    double Bx{ 0 }, By{ 0 }, Bz{ 0 }, B{ 0 };
    double ElementForcex{ 0 }, ElementForcey{ 0 }, ElementForcez{ 0 }, ElementForce{ 0 };
    double RHSContri{ 0 };
    bool linearflag{ true };
};

int m_num_nodes;
int m_num_vtxele;
int m_num_edgele;
int m_num_triele;
int num_freenodes;
CNode* mp_node;
CVtxElement* mp_vtxele;
CEdgElement* mp_edgele;
CTriElement* mp_triele;

CNode* d_mp_node;
CTriElement* d_mp_triele;

void FEM_Host_Data_Prepare();
void LoadMeshInfo();
void processBoundaryCondition();
void processLoad();
void processNDDRNode();
void makeTrangle();
void processMaterial();
void GPUInitialMallocCopy();
void writeVtkFile();
__global__ void UpdateSolutiontoA1(int m_num_nodes, CTriElement* d_MyElem, CNode* d_MyNode);
__global__ void UpdateAttoAtold(int m_num_nodes, CNode* d_MyNode);
__global__ void calculateGlobalError(int m_num_nodes, CTriElement* d_MyElem, CNode* d_MyNode, double* a, double* b, double* d_error);

__device__ double getV(double B);
__device__ double getdVdB2(double B);

void Free();
int main()
{
    FEM_Host_Data_Prepare();
    GPUInitialMallocCopy();

    cudaEvent_t start, stop;//unit: ms

    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    double* a, * b, error, * d_error;
    cudaMalloc(&a, m_num_nodes * sizeof(double));
    cudaMalloc(&b, m_num_nodes * sizeof(double));
    cudaMalloc(&d_error, sizeof(double));
    for (int iter = 0; iter < 671; ++iter) {
        UpdateSolutiontoA1 << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_triele, d_mp_node);
        //calculateGlobalError << <256, 256 >> > (m_num_nodes, d_mp_triele, d_mp_node, a, b, d_error);
        //cudaDeviceSynchronize();
        //cudaMemcpy(&error, d_error, sizeof(double), cudaMemcpyDeviceToHost);
        //printf("iter: %d, error: %.20f\n",iter, error);
        UpdateAttoAtold << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_node);
    }
    cudaDeviceSynchronize();
    //GET_TIME(T_end);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    //printf("\nElapsed time for NDDR Iteration is %.4f s\n", (T_end - T_start));
    printf("\nElapsed time for NDDR Iteration is %.4f s\n", elapsedTime / 1000);
    cudaMemcpy(mp_node, d_mp_node, m_num_nodes * sizeof(CNode), cudaMemcpyDeviceToHost);
    cudaMemcpy(mp_triele, d_mp_triele, m_num_triele * sizeof(CTriElement), cudaMemcpyDeviceToHost);
    writeVtkFile();

    Free();
    return 0;
}


void FEM_Host_Data_Prepare() {
    LoadMeshInfo();
    processBoundaryCondition();
    processLoad();
    processNDDRNode();
    makeTrangle();
    processMaterial();


}

void LoadMeshInfo() {
    string meshfile = "../model/model1848.mphtxt";
    //分网文件为mphtxt的情况
    
    cout << meshfile << endl;
    FILE* fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, meshfile.c_str(), "r");
    if (!fp) {
        std::cout << "Error: openning file!" << endl;
        exit(0);
    }

    //--------------read the head-----------------------------
    for (int i = 0; i < 17; i++) {
        fgets(ch, 256, fp);
    }
    //-----------------mesh point-----------------------------
    if (fscanf_s(fp, "%d # number of mesh vertices\n", &m_num_nodes) != 1) {
        std::cout << "error: reading num_bdr_ns!" << endl;
        exit(0);
    }
    else std::cout << m_num_nodes << "number of mesh points." << endl;
    mp_node = new CNode[m_num_nodes];
    int pts_ind;//the beginning of the points index

    if (fscanf_s(fp, "%d # lowest mesh vertex index\n", &pts_ind) != 1) {
        std::cout << "error: reading pts_ind!" << endl;
        exit(0);
    }
    fgets(ch, 256, fp);

    for (int i = pts_ind; i < m_num_nodes; i++) {
        //读取x,y坐标
        if (fscanf_s(fp, "%lf %lf \n", &(mp_node[i].x), &(mp_node[i].y)) != 2) {
            std::cout << "error: reading mesh point!" << endl;
            exit(0);
        }
        //else{
        //    std::cout << mp_2DNode[i].x << " " << mp_2DNode[i].y << "\n";
        //}
    }
    //---------------vertexnode-------------------------------
    for (int i = 0; i < 7; i++)
        fgets(ch, 256, fp);
    int num_vtx_ns;

    if (fscanf_s(fp, "%d # number of vertices per element\n", &num_vtx_ns) != 1) {
        std::cout << "error: reading num_vtx_ns!" << endl;
        exit(0);
    }

    if (fscanf_s(fp, "%d # number of elements\n", &m_num_vtxele) != 1) {
        std::cout << "error: reading m_num_vtxele!" << endl;
        exit(0);
    }
    else std::cout << m_num_vtxele << "number of vertex elements." << endl;
    fgets(ch, 256, fp);
    //    std::cout << m_num_vtxele;
    mp_vtxele = new CVtxElement[m_num_vtxele];
    for (int i = 0; i < m_num_vtxele; i++) {
        if (fscanf_s(fp, "%d \n", &((mp_vtxele + i)->n)) != 1) {
            std::cout << "error: reading vertex element points!" << endl;
            exit(0);
        }
    }
    //---------------vertexdomain-------------------------------
    for (int i = 0; i < 2; i++)
        fgets(ch, 256, fp);
    for (int i = 0; i < m_num_vtxele; i++) {
        if (fscanf_s(fp, "%d \n", &(mp_vtxele[i].domain)) != 1) {
            std::cout << "error: reading vertex domain!" << endl;
            exit(0);
        }
        else {
            mp_vtxele[i].domain++;
        }
    }
    //----------------edgnode-----------------------------------
    for (int i = 0; i < 6; i++) {
        fgets(ch, 256, fp);
    }
    if (fscanf_s(fp, "%d # number of elements\n", &m_num_edgele) != 1) {
        std::cout << "error: reading m_num_edgele" << endl;
        exit(0);
    }
    else std::cout << m_num_edgele << "number of edge elements." << endl;
    mp_edgele = new CEdgElement[m_num_edgele];
    fgets(ch, 256, fp);
    for (int i = 0; i < m_num_edgele; i++) {
        if (fscanf_s(fp, "%d %d\n", &mp_edgele[i].n[0], &mp_edgele[i].n[1]) != 2) {
            std::cout << "error: reading edg element points!" << endl;
            exit(0);
        }
    }
    //----------------edgdomain---------------------------------
    for (int i = 0; i < 2; i++) {
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_edgele; i++) {
        if (fscanf_s(fp, "%d \n", &mp_edgele[i].domain) != 1) {
            std::cout << "error: reading edgdomain!" << endl;
            exit(0);
        }
        else {
            mp_edgele[i].domain++;
            //            std::cout << "mp_edgele: " << mp_edgele[i].domain;
        }
    }
    //----------------trinode-----------------------------------
    for (int i = 0; i < 6; i++) {
        fgets(ch, 256, fp);
    }
    if (fscanf_s(fp, "%d # number of elements\n", &m_num_triele) != 1) {
        std::cout << "error: reading m_num_triele!" << endl;
        exit(0);
    }
    else std::cout << m_num_triele << "number of triangle elements." << endl;
    fgets(ch, 256, fp);
    mp_triele = new CTriElement[m_num_triele];
    for (int i = 0; i < m_num_triele; i++) {
        if (fscanf_s(fp, "%d %d %d \n", &mp_triele[i].n[0], &mp_triele[i].n[1], &mp_triele[i].n[2]) != 3) {
            std::cout << "error: reading elements points!" << endl;
            exit(0);
        }
    }
    //----------------tridomain---------------------------------
    for (int i = 0; i < 2; i++) {
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_triele; i++) {
        if (fscanf_s(fp, "%d \n", &mp_triele[i].domain) != 1) {
            std::cout << "error: reading tridomain!" << endl;
            exit(0);
        }
        else {
            //mp_triele[i].domain++;
            //cout << mp_triele[i].domain << endl;
        }
    }
    fclose(fp);

    std::cout << "load2dmeshcomsol successfully!" << endl;;
}

void processBoundaryCondition() {
    //这里的boundarymap只用了第一类边界条件，只考虑线单元域
//1.提取线所在全部域，2.提取域对应的单元编号，3.提取节点
//4.去重 5.给node设置边界条件 6.将边界点排序到末尾
    num_freenodes = m_num_nodes;
    std::vector<int> boundarynodes;
    std::vector<int> boundarylist{ 1, 2, 3, 5, 7, 31, 32, 37, 38, 48, 49, 52, 53, 59, 60, 61, 62, 63, 64 };
 
    //原有版本的边界检索
    for (auto a : boundarylist) {
        int boundarynodedomain = a;
        for (int i = 0; i < m_num_edgele; ++i) {
            if (mp_edgele[i].domain == boundarynodedomain) {
                boundarynodes.push_back(mp_edgele[i].n[0]);
                boundarynodes.push_back(mp_edgele[i].n[1]);
            }
        }
    }

    //将x=0的全部节点标记为边界节点
    for (int i = 0; i < m_num_nodes; ++i) {
        if (mp_node[i].x == 0) {
            boundarynodes.push_back(i);
        }
    }

    std::sort(boundarynodes.begin(), boundarynodes.end());
    boundarynodes.erase(unique(boundarynodes.begin(), boundarynodes.end()), boundarynodes.end());

    ////输出边界点，可以通过matlab绘图判断是否检索出全部的边界
    //for (auto a :boundarynodes) {
    //	cout << mp_node[a].x << " " << mp_node[a].y << endl;
    //}

    num_freenodes -= boundarynodes.size();
    for (auto a : boundarynodes) {
        mp_node[a].bdr = 1;
    }
}

void processLoad() {
    //设置单元负载
    for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
        int domain = mp_triele[i_tri].domain;
        if (domain == 5) {
            mp_triele[i_tri].J = CurrentDensity;
        }
    }
}

void processNDDRNode()
{
    for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
        for (int j = 0; j < 3; ++j) {
            int n = mp_triele[i_tri].n[j];
            int id = mp_node[n].NumberofNeighbourElement;
            mp_node[n].NeighbourElementId[id] = i_tri;
            mp_node[n].NeighbourElementNumber[id] = j;
            mp_node[n].NumberofNeighbourElement++;
            //printf("global n:%d, NumberofNeighbourElement:%d\n", n, mp_node[n].NumberofNeighbourElement);
        }
        //printf("global ele:%d, J:%f\n", i_tri, mp_triele[i_tri].J);
    }
}

void makeTrangle() {
    for (int index = 0; index < m_num_triele; ++index) {
        int k, m, n;
        double q0, q1, q2, r0, r1, r2, area;
        k = mp_triele[index].n[0];
        m = mp_triele[index].n[1];
        n = mp_triele[index].n[2];

        q0 = mp_node[m].y - mp_node[n].y;
        mp_triele[index].Q[0] = q0;
        q1 = mp_node[n].y - mp_node[k].y;
        mp_triele[index].Q[1] = q1;
        q2 = mp_node[k].y - mp_node[m].y;
        mp_triele[index].Q[2] = q2;

        r0 = mp_node[n].x - mp_node[m].x;
        mp_triele[index].R[0] = r0;
        r1 = mp_node[k].x - mp_node[n].x;
        mp_triele[index].R[1] = r1;
        r2 = mp_node[m].x - mp_node[k].x;
        mp_triele[index].R[2] = r2;

        area = 0.5 * std::abs(q1 * r2 - r1 * q2);
        mp_triele[index].area = area;

        mp_triele[index].rc = (mp_node[k].x +
            mp_node[m].x +
            mp_node[n].x) / 3;
        mp_triele[index].zc = (mp_node[k].y +
            mp_node[m].y +
            mp_node[n].y) / 3;

        int flag = 0;
        for (int f = 0; f < 3; f++) {
            if (mp_node[mp_triele[index].n[f]].x < 1e-7) {
                flag++;
            }
        }

        //计算三角形重心半径
        if (flag == 2) {
            mp_triele[index].xdot = mp_triele[index].rc;
        }
        else {
            mp_triele[index].xdot = 1 / (mp_node[k].x + mp_node[m].x);
            mp_triele[index].xdot += 1 / (mp_node[k].x + mp_node[n].x);
            mp_triele[index].xdot += 1 / (mp_node[m].x + mp_node[n].x);
            mp_triele[index].xdot = 1.5 / mp_triele[index].xdot;
        }

        //计算一阶三角形轴对称单元系数矩阵
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                mp_triele[index].C[i][j] = ((mp_triele[index].R[i] * mp_triele[index].R[j] + mp_triele[index].Q[i] * mp_triele[index].Q[j])) / (4.0 * mp_triele[index].area);
            }
        }
    }
}

void processMaterial() {
    for (int i = 0; i < m_num_triele; ++i) {
        if (mp_triele[i].domain == 3 || mp_triele[i].domain == 4) {
            mp_triele[i].linearflag = false;
        }
        else {
            mp_triele[i].linearflag = true;
        }
    }
}

void GPUInitialMallocCopy() {
    //cuda initiallze
    int num_devices, device;
    cudaGetDeviceCount(&num_devices);
    if (num_devices > 1) {
        int max_multiprocessors = 0, max_device = 0;
        for (device = 0; device < num_devices; device++) {
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, device);
            if (max_multiprocessors < properties.multiProcessorCount) {
                max_multiprocessors = properties.multiProcessorCount;
                max_device = device;
            }
        }

    }
    cudaSetDevice(0);

    //分配节点和TriElement的内存
    cudaMalloc(&d_mp_node, m_num_nodes * sizeof(CNode));
    cudaMemcpy(d_mp_node, mp_node, m_num_nodes * sizeof(CNode), cudaMemcpyHostToDevice);
    cudaMalloc(&d_mp_triele, m_num_triele * sizeof(CTriElement));
    cudaMemcpy(d_mp_triele, mp_triele, m_num_triele * sizeof(CTriElement), cudaMemcpyHostToDevice);
}

void Free() {
    cudaFree(d_mp_node);
    cudaFree(d_mp_triele);
    delete[] mp_node;
    delete[] mp_vtxele;
    delete[] mp_edgele;
    delete[] mp_triele;
}

void writeVtkFile()
{
    //处理数据
    for (int i = 0; i < m_num_nodes; ++i) {
        if (mp_node[i].x != 0)
            mp_node[i].A = mp_node[i].At / mp_node[i].x;
    }
    char name[] = "result.vtk";
    FILE* fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, name, "w");
    if (!fp) {
        std::cout << "Error: openning file!" << endl;
        exit(0);
    }
    /*
         1: points
         3: line
         5: Triangular element
         9: Quadrilateral element
        10: Tetrahedral element
        12: Hexahedral element
        13: Triangular prism element
        14: Pyramid element
    */
    /** 数据版本声明 **/
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    /** 标题 **/
    fprintf(fp, "vtk title\n");
    /** 文件格式声明 **/
    fprintf(fp, "ASCII\n");
    /** 几何拓扑结构 **/
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    //节点
    fprintf(fp, "\nPOINTS %d float\n", m_num_nodes);
    for (int i = 0; i < m_num_nodes; ++i) {
        fprintf(fp, "%lf %lf %lf\n", mp_node[i].x, mp_node[i].y, mp_node[i].z);
    }
    //一阶三角形单元
    fprintf(fp, "\nCELLS %d %d\n", m_num_triele, 4 * m_num_triele);
    for (int i = 0; i < m_num_triele; ++i) {
        fprintf(fp, "3 %d %d %d\n", mp_triele[i].n[0], mp_triele[i].n[1], mp_triele[i].n[2]);
    }
    fprintf(fp, "\nCELL_TYPES %d\n", m_num_triele);
    int type = 5;
    for (int i = 0; i < m_num_triele; ++i) {
        fprintf(fp, "%d\n", type);
    }
    //节点磁矢位
    fprintf(fp, "\nPOINT_DATA %d\n", m_num_nodes);
    fprintf(fp, "SCALARS A double 1\n");
    fprintf(fp, "LOOKUP_TABLE %s\n", "Atable");
    for (int i = 0; i < m_num_nodes; ++i) {
        fprintf(fp, "%f\n", mp_node[i].A);
    }

    ////单元标量磁感应强度
    //fprintf(fp, "\nCELL_DATA %d\n", m_num_triele);
    //fprintf(fp, "SCALARS %s double %d\n", "Bnorm", 1);
    //fprintf(fp, "LOOKUP_TABLE %s\n", "Btable");
    //for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
    //    fprintf(fp, "%lf\n", mp_triele[i_tri].B);
    //}

    //fprintf(fp, "\nVECTORS %s double\n", "Bvector");
    //for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
    //    fprintf(fp, "%lf %lf %lf\n", mp_triele[i_tri].Bx, mp_triele[i_tri].By, 0.0);
    //}
    //fclose(fp);
}

__global__ void UpdateSolutiontoA1(int m_num_nodes, CTriElement* d_MyElem, CNode* d_MyNode) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= m_num_nodes)
        return;

    double A, At, dAt, Jac, k, J0, k0, NRerror, Js;
    double B2, B, V, Vt, VB2, B2A;
    int ID;
    int RelaxCount, count = 0;
    double RelaxFactor = 1;
    double AtLocal[3]{ 0 ,0, 0 }, ALocal[3]{ 0, 0, 0 };
    int NRCount = 6;
    if (d_MyNode[i].bdr != 1)
    {
        At = 0; dAt = 0; NRerror = 1; count = 0;
        while (count < NRCount)
        {
            Jac = 0; k = 0;
            for (int j = 0; j < d_MyNode[i].NumberofNeighbourElement; j++)
            {
                ID = d_MyNode[i].NeighbourElementId[j];
                Js = d_MyElem[ID].J * d_MyElem[ID].area / 3;
                int nodenumber = d_MyNode[i].NeighbourElementNumber[j];
                for (int m = 0; m < 3; ++m) {
                    if (m == nodenumber) {
                        AtLocal[m] = At;
                    }
                    else {
                        AtLocal[m] = d_MyNode[d_MyElem[ID].n[m]].At_old;
                    }
                    ALocal[m] = AtLocal[m] / d_MyElem[ID].xdot;
                    A = At / d_MyElem[ID].xdot;
                }
                double RHSContri = 0;
                for (int m = 0; m < 3; ++m) {

                    RHSContri += d_MyElem[ID].C[nodenumber][m] * AtLocal[m];
                }

                //线性单元，计算右侧列向量时，也要考虑A？？？
                if (d_MyElem[ID].linearflag == true)	//线性空气单元
                {
                    Vt = 1.0 / (4 * 3.14159 * 1e-7) / d_MyElem[ID].xdot;
                    J0 = Vt * d_MyElem[ID].C[nodenumber][nodenumber];
                    k0 = Js - Vt * RHSContri;
                    //V = 1.0 / d_MyElem[ID].material->getMu(0);
                    //J0 = V * d_MyElem[ID].C[nodenumber][nodenumber];
                    //k0 = Js - V * RHSContri;
                }

                else	//非线性单元
                {
                    B2 = -1 / d_MyElem[ID].area * (d_MyElem[ID].C[0][1] * (ALocal[0] - ALocal[1]) * (ALocal[0] - ALocal[1]) + d_MyElem[ID].C[1][2] * (ALocal[1] - ALocal[2]) * (ALocal[1] - ALocal[2]) + d_MyElem[ID].C[0][2] * (ALocal[0] - ALocal[2]) * (ALocal[0] - ALocal[2]));
                    B = sqrt(B2);	//计算B，用于处理非线性
                    //Vt = getV(B) / d_MyElem[ID].xdot;
                    //VB2 = getdVdB2(B);
                    if (B <= 0.6) {
                        Vt = 500;
                    }
                    else {
                        Vt =  500 + 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
                    }
                    Vt = Vt / d_MyElem[ID].xdot;
                    if (B <= 0.6) {
                        VB2 =  0;
                    }
                    else {
                        VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
                    }
                    B2A = -1 / d_MyElem[ID].area * (2 * d_MyElem[ID].C[nodenumber][0] * (A - ALocal[0]) + 2 * d_MyElem[ID].C[nodenumber][1] * (A - ALocal[1]) + 2 * d_MyElem[ID].C[nodenumber][2] * (A - ALocal[2]));
                    J0 = B2A * VB2 * RHSContri / d_MyElem[ID].xdot / d_MyElem[ID].xdot + Vt * d_MyElem[ID].C[nodenumber][nodenumber];
                    k0 = Js - Vt * RHSContri;

                    //B2 = -1 / d_MyElem[ID].area * (d_MyElem[ID].C[0][1] * (AtLocal[0] - AtLocal[1]) * (AtLocal[0] - AtLocal[1]) + d_MyElem[ID].C[1][2] * (AtLocal[1] - AtLocal[2]) * (AtLocal[1] - AtLocal[2]) + d_MyElem[ID].C[0][2] * (AtLocal[0] - AtLocal[2]) * (AtLocal[0] - AtLocal[2]));
                    //B = sqrt(B2);	//计算B，用于处理非线性
                    //V = 1.0 / d_MyElem[ID].material->getMu(B);
                    //VB2 = d_MyElem[ID].material->getdvdB2(B);
                    //B2A = -1 / d_MyElem[ID].area * (2 * d_MyElem[ID].C[nodenumber][0] * (A - AtLocal[0]) + 2 * d_MyElem[ID].C[nodenumber][1] * (A - AtLocal[1]) + 2 * d_MyElem[ID].C[nodenumber][2] * (A - AtLocal[2]));
                    //J0 = B2A * VB2 * RHSContri + V * d_MyElem[ID].C[nodenumber][nodenumber];
                    //k0 = Js - V * RHSContri;
                }

                Jac = Jac + J0;
                k = k + k0;
            }
            dAt = k / Jac;
            At = At + RelaxFactor * dAt;

            double a = dAt * dAt;
            double b = At * At;
            double NRerror = sqrt(a) / sqrt(b);
            if (At == 0) {
                break;
            }
            if (NRerror > 1e-9) {
                count++;
                continue;
            }
            else {
                //if (NRiter != 1) {
                //	printf("NRerror: %.20f\n", NRerror);
                //	cout << "NRerror: " << NRerror << endl;
                //	printf("n: %d, NRiter: %d\n", n, NRiter);
                //}
                break;
            }
            //count++;
        }

        d_MyNode[i].At = At;
    }

    ////判断全局收敛性
    //double error = 0, a = 0, b = 0;
    //for (int i = 0; i < m_num_nodes; ++i) {
    //    a += (d_MyNode[i].At - d_MyNode[i].At_old) * (d_MyNode[i].At - d_MyNode[i].At_old);
    //    b += d_MyNode[i].At * d_MyNode[i].At;
    //}
    //error = sqrt(a) / sqrt(b);
    ////if ((iter + 1) % 100 == 0) {
    ////	cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
    ////}
    //cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
    //if (error > maxerror) {
    //    //这一步影响效率
    //    for (int i = 0; i < m_num_nodes; ++i) {
    //        d_MyNode[i].At_old = d_MyNode[i].At;
    //    }
    //}
    //else {
    //    cout << "Iteration step: " << iter + 1 << endl;
    //    cout << "Nonlinear NDDR iteration finish.\n";
    //    return;
    //}

    //for (int i = 0; i < m_num_nodes; ++i) {
    //    if (d_MyNode[i].x != 0)
    //        d_MyNode[i].A = d_MyNode[i].At / d_MyNode[i].x;
    //}
}


__global__ void UpdateAttoAtold(int m_num_nodes, CNode* d_MyNode) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= m_num_nodes)
        return;
    d_MyNode[i].At_old = d_MyNode[i].At;
}

__global__ void calculateGlobalError(int m_num_nodes, CTriElement* d_MyElem, CNode* d_MyNode, double* a, double* b, double* d_error) {
    int n = blockDim.x * blockIdx.x + threadIdx.x;
    if (n >= m_num_nodes) {
        return;
    }

    a[n] = (d_MyNode[n].At - d_MyNode[n].At_old) * (d_MyNode[n].At - d_MyNode[n].At_old);
    b[n] = d_MyNode[n].At - d_MyNode[n].At_old;
    //d_error = 0;
    __syncthreads();

    //a和b归约求和
    int leng = m_num_nodes;
    for (int i = m_num_nodes / 2.0 + 0.5; i > 1; i = i / 2.0 + 0.5) {
        if (n < i)
        {

            if (n + i < leng)
            {
                a[n] += a[n + i];
                b[n] += b[n + i];
            }
        }
        __syncthreads();
        leng = leng / 2.0 + 0.5;
    }

    if (n == 0) {
        a[0] = a[0] + a[1];
        b[0] = b[0] + b[1];
        *d_error = sqrtf(a[0]) / sqrtf(b[0]);
    }
}

__device__ double getV(double B) {
    if (B <= 0.6) {
        return 500;
    }
    else {
        return 500 + 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
    }
}

__device__ double getdVdB2(double B) {
    if (B <= 0.6) {
        return 0;
    }
    else {
        return (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
    }
}