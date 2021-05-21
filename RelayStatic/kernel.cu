#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <Algorithm>
#include <string>
#define Mu0 1.256637e-6;
#define PI 3.1415927;

const int CudaThrdNum = 128;
const int CudaBlckNum = 128;
const double CurrentDensity = 8381893.016;
const int maxitersteps = 20000;

const double maxerror = 1e-9;
__device__ double Bdata[] = { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4 };
__device__ double Hdata[] = { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
__device__ int BHpoints = 16;
using namespace std;

//运动信息
double current[101]/* = { 0.000000, 0.017863, 0.035336, 0.052028, 0.068390, 0.084379, 0.099983, 0.116192, 0.128333, 0.143884, 0.158808, 0.170294, 0.185229, 0.197768, 0.209873, 0.221648, 0.233348, 0.245268, 0.257065, 0.268852, 0.280279, 0.291660, 0.302664, 0.313246, 0.323835, 0.334095, 0.344348, 0.354363, 0.364181, 0.373710, 0.381755, 0.390556, 0.398050, 0.406149, 0.412667, 0.419986, 0.423393, 0.426914, 0.427973, 0.427123, 0.419883, 0.420140, 0.380406, 0.394262, 0.322280, 0.264271, 0.275031, 0.292164, 0.309064, 0.326167, 0.342338, 0.357881, 0.372799, 0.387087, 0.400791, 0.413996, 0.426584, 0.438703, 0.449832, 0.460407, 0.470367, 0.479662, 0.488457, 0.496547, 0.504091, 0.510959, 0.517875, 0.523941, 0.529457, 0.534735, 0.539611, 0.544017, 0.548248, 0.552263, 0.556181, 0.559601, 0.562426, 0.565602, 0.568039, 0.570573, 0.573147, 0.575426, 0.577307, 0.578795, 0.580537, 0.582063, 0.583529, 0.585076, 0.585960, 0.587053, 0.588023, 0.589394, 0.590167, 0.590551, 0.591922, 0.592885, 0.593365, 0.593414, 0.594172, 0.594799, 0.594968 }*/;
double dis[101];
double velocity[101];
double acc[101];
double magneticforce[101];
double springforce[101];
double flux[101]/* = { 0.000000, 0.011649, 0.022951, 0.033938, 0.044603, 0.054956, 0.065021, 0.074693, 0.084189, 0.093420, 0.102245, 0.110904, 0.119200, 0.127245, 0.135048, 0.142616, 0.149951, 0.157046, 0.163906, 0.170529, 0.176928, 0.183095, 0.189041, 0.194784, 0.200307, 0.205633, 0.210754, 0.215677, 0.220404, 0.224941, 0.229328, 0.233541, 0.237686, 0.241607, 0.245419, 0.249023, 0.252580, 0.256066, 0.259562, 0.263109, 0.266877, 0.270605, 0.274998, 0.279217, 0.284745, 0.291486, 0.297971, 0.304115, 0.309954, 0.315417, 0.320563, 0.325403, 0.329946, 0.334205, 0.338211, 0.341956, 0.345449, 0.348677, 0.351687, 0.354486, 0.357079, 0.359496, 0.361727, 0.363805, 0.365736, 0.367543, 0.369194, 0.370715, 0.372143, 0.373465, 0.374692, 0.375837, 0.376902, 0.377884, 0.378781, 0.379608, 0.380387, 0.381094, 0.381762, 0.382382, 0.382939, 0.383431, 0.383885, 0.384330, 0.384742, 0.385126, 0.385476, 0.385775, 0.386080, 0.386378, 0.386648, 0.386882, 0.387105, 0.387347, 0.387533, 0.387675, 0.387833, 0.387995, 0.388140, 0.388273, 0.388423}*/;
double mass = 0.024;
double h = 5e-4;
double U = 24;
double R = 40;
double* d_FluxArray;
double* d_Flux;
double position[] = { 0, 0.0016999, 0.0017, 0.0027 };
double force[] = { -6.0, -6.63, -13.63, -27.0 };
int staticsteps;
int dynamicsteps;

struct CNode
{
    double x{ 0 }, y{ 0 }, z{ 0 };
    int bdr{ 0 }; //�߽�����
    double A{ 0 };
    double A_old{ 0 };
    double At{ 0 };
    double At_old{ 0 };
    int NumberofNeighbourElement{ 0 };
    int NeighbourElementId[15];	//�ͽڵ���صĵ�Ԫ���
    int NeighbourElementNumber[15];	//�ڵ��ڶ�Ӧ��Ԫ�еı��
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
    double C[3][3];// ��Ԫϵ������
    double area{ 0 }; 
    double rc, zc;
    double xdot;
    int domain{ 0 };
    double J{ 0 };  //���أ���ʱֻ���ǵ�����ֱ����������Ҫ������һ����
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
void LoadMeshInfo(string meshfile);
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
__global__ void updateLoad(int m_num_triele, CTriElement* d_MyElem, double current);
__global__ void solveFlux(int m_num_triele, CTriElement* d_MyElem, CNode* d_MyNode, double* d_FluxArray, double* d_Flux);

__device__ double getV(double B);
__device__ double getdVdB2(double B);
__device__ double getkHb(double B, double* k, double* H, double* b);

void Free();

void solvestatic();
void solvedynamic();
void solveWeakCouple(int step);
double solveSpringForce(double pos);
double solveMagneticForce();

int main()
{
    //solvestatic();
    solvedynamic();
    return 0;
}


void FEM_Host_Data_Prepare() {
    LoadMeshInfo("D:/femplatform/model/geo/modelcomsol_dynamic_NDDRGPU/modelwithband_0.mphtxt");
    processBoundaryCondition();
    processLoad();
    processNDDRNode();
    makeTrangle();
    processMaterial();


}

void LoadMeshInfo(string meshfile) {
    //string meshfile = "D:/femplatform/model/geo/modelcomsol_dynamic_NDDRGPU/modelwithband_0.mphtxt";
    //string meshfile = "D:/femplatform/model/model1848.mphtxt";
    //�����ļ�Ϊmphtxt�����
    
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
        //��ȡx,y����
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
    //�����boundarymapֻ���˵�һ��߽�������ֻ�����ߵ�Ԫ��
//1.��ȡ������ȫ����2.��ȡ���Ӧ�ĵ�Ԫ��ţ�3.��ȡ�ڵ�
//4.ȥ�� 5.��node���ñ߽����� 6.���߽������ĩβ
    num_freenodes = m_num_nodes;
    std::vector<int> boundarynodes;
    std::vector<int> boundarylist{ 1, 2, 3, 5, 6, 8, 40, 43 };
 
    //ԭ�а汾�ı߽����
    for (auto a : boundarylist) {
        int boundarynodedomain = a;
        for (int i = 0; i < m_num_edgele; ++i) {
            if (mp_edgele[i].domain == boundarynodedomain) {
                boundarynodes.push_back(mp_edgele[i].n[0]);
                boundarynodes.push_back(mp_edgele[i].n[1]);
            }
        }
    }

    //��x=0��ȫ���ڵ���Ϊ�߽�ڵ�
    for (int i = 0; i < m_num_nodes; ++i) {
        if (mp_node[i].x == 0) {
            boundarynodes.push_back(i);
        }
    }

    std::sort(boundarynodes.begin(), boundarynodes.end());
    boundarynodes.erase(unique(boundarynodes.begin(), boundarynodes.end()), boundarynodes.end());

    ////����߽�㣬����ͨ��matlab��ͼ�ж��Ƿ������ȫ���ı߽�
    //for (auto a :boundarynodes) {
    //	cout << mp_node[a].x << " " << mp_node[a].y << endl;
    //}

    num_freenodes -= boundarynodes.size();
    for (auto a : boundarynodes) {
        mp_node[a].bdr = 1;
    }
}

void processLoad() {
    //���õ�Ԫ����
    for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
        int domain = mp_triele[i_tri].domain;
        if (domain == 7) {
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

        //�������������İ뾶
        if (flag == 2) {
            mp_triele[index].xdot = mp_triele[index].rc;
        }
        else {
            mp_triele[index].xdot = 1 / (mp_node[k].x + mp_node[m].x);
            mp_triele[index].xdot += 1 / (mp_node[k].x + mp_node[n].x);
            mp_triele[index].xdot += 1 / (mp_node[m].x + mp_node[n].x);
            mp_triele[index].xdot = 1.5 / mp_triele[index].xdot;
        }

        //����һ����������ԳƵ�Ԫϵ������
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                mp_triele[index].C[i][j] = ((mp_triele[index].R[i] * mp_triele[index].R[j] + mp_triele[index].Q[i] * mp_triele[index].Q[j])) / (4.0 * mp_triele[index].area);
            }
        }
    }
}

void processMaterial() {
    for (int i = 0; i < m_num_triele; ++i) {
        if (mp_triele[i].domain == 4 || mp_triele[i].domain == 5) {
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

    //����ڵ��TriElement���ڴ�
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
    //��������
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
    /** ���ݰ汾���� **/
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    /** ���� **/
    fprintf(fp, "vtk title\n");
    /** �ļ���ʽ���� **/
    fprintf(fp, "ASCII\n");
    /** �������˽ṹ **/
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    //�ڵ�
    fprintf(fp, "\nPOINTS %d float\n", m_num_nodes);
    for (int i = 0; i < m_num_nodes; ++i) {
        fprintf(fp, "%lf %lf %lf\n", mp_node[i].x, mp_node[i].y, mp_node[i].z);
    }
    //һ�������ε�Ԫ
    fprintf(fp, "\nCELLS %d %d\n", m_num_triele, 4 * m_num_triele);
    for (int i = 0; i < m_num_triele; ++i) {
        fprintf(fp, "3 %d %d %d\n", mp_triele[i].n[0], mp_triele[i].n[1], mp_triele[i].n[2]);
    }
    fprintf(fp, "\nCELL_TYPES %d\n", m_num_triele);
    int type = 5;
    for (int i = 0; i < m_num_triele; ++i) {
        fprintf(fp, "%d\n", type);
    }
    //�ڵ��ʸλ
    fprintf(fp, "\nPOINT_DATA %d\n", m_num_nodes);
    fprintf(fp, "SCALARS A double 1\n");
    fprintf(fp, "LOOKUP_TABLE %s\n", "Atable");
    for (int i = 0; i < m_num_nodes; ++i) {
        fprintf(fp, "%f\n", mp_node[i].A);
    }

    ////��Ԫ�����Ÿ�Ӧǿ��
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

void solvestatic() {
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
    for (int iter = 0; iter < 30000; ++iter) {
        UpdateSolutiontoA1 << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_triele, d_mp_node);
        calculateGlobalError << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_triele, d_mp_node, a, b, d_error);
        cudaDeviceSynchronize();
        cudaMemcpy(&error, d_error, sizeof(double), cudaMemcpyDeviceToHost);
        printf("iter: %d, error: %.20f\n", iter, error);
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
}

void solvedynamic() {
    FEM_Host_Data_Prepare();
    GPUInitialMallocCopy();

    cudaMalloc(&d_FluxArray, m_num_triele * sizeof(double));
    cudaMalloc(&d_Flux, sizeof(double));

    dis[0] = 0;
    velocity[0] = 0;
    acc[0] = 0;
    springforce[0] = solveSpringForce(0);
    magneticforce[0] = 0;
    current[0] = 0;
    flux[0] = 0;
    int dynamicstepsarray[101];
    double time[101];
    dynamicstepsarray[0] = 0;
    time[0] = 0;


    bool stopflag = false;
    clock_t start, end;
    start = clock();
    for (int i = 1; i < 101; ++i) {
        dynamicsteps = 0;
        clock_t start, end;
        start = clock();
        cout << "solve step " << i << "...\n";
        acc[i] = (magneticforce[i - 1] + springforce[i - 1]) / mass;
        if (acc[i] < 0) acc[i] = 0;
        velocity[i] = velocity[i - 1] + h * acc[i - 1];
        dis[i] = dis[i - 1] + h * velocity[i - 1];

        if (dis[i] >= 0.0024) {
            dis[i] = 0.0024;
            if (stopflag == false) {
                stopflag = true;
            }
            else {
                velocity[i] = 0;
                acc[i] = 0;
            }
        }

        //后向欧拉法计算
        springforce[i] = solveSpringForce(dis[i]);
        acc[i] = (magneticforce[i - 1] + springforce[i - 1]) / mass;
        if (acc[i] < 0) acc[i] = 0;
        velocity[i] = velocity[i - 1] + h * acc[i];
        dis[i] = dis[i - 1] + h * velocity[i];

        if (dis[i] >= 0.0024) {
            dis[i] = 0.0024;
            if (stopflag == false) {
                stopflag = true;
            }
            else {
                velocity[i] = 0;
                acc[i] = 0;
            }
        }

        if (i >= 30 && i <= 46) {
            string meshfile = "D:/femplatform/model/geo/modelcomsol_dynamic_NDDRGPU/modelwithband_";
            meshfile += to_string(i) + ".mphtxt";
            Free();
            LoadMeshInfo(meshfile);
            processBoundaryCondition();
            processLoad();
            processNDDRNode();
            makeTrangle();
            processMaterial();
            cudaMalloc(&d_mp_node, m_num_nodes * sizeof(CNode));
            cudaMemcpy(d_mp_node, mp_node, m_num_nodes * sizeof(CNode), cudaMemcpyHostToDevice);
            cudaMalloc(&d_mp_triele, m_num_triele * sizeof(CTriElement));
            cudaMemcpy(d_mp_triele, mp_triele, m_num_triele * sizeof(CTriElement), cudaMemcpyHostToDevice);

        }
        solveWeakCouple(i);
        magneticforce[i] = solveMagneticForce();
        springforce[i] = solveSpringForce(dis[i]);
        dynamicstepsarray[i] = dynamicsteps;
        end = clock();
        time[i] = double(end - start) / CLOCKS_PER_SEC;
        printf("dynamicsteps: %d, time: %f\n", dynamicstepsarray[i], time[i]);
        printf("step: %d, dis: %f, velocity: %f, acc: %f, springforce: %f, magneticforce: %f\n\n", i, dis[i], velocity[i], acc[i], springforce[i], magneticforce[i]);
    }

    //写入结果文件
    char fn[256];
    sprintf(fn, "%s.m", "RelayDynamic");
    FILE* fp;
    fp = fopen(fn, "w");
    fprintf(fp, "%%output by FEEM\n");
    fprintf(fp, "%%timesteps, displacements, velocities, accelerations, magneticforce, current, flux, steps, time\n");

    fprintf(fp, "results = [\n");
    for (int i = 0; i < 101; ++i) {
        fprintf(fp, "%10.8e,", i * h);
        fprintf(fp, "%10.8e,", dis[i]);
        fprintf(fp, "%10.8e,", velocity[i]);
        fprintf(fp, "%10.8e,", acc[i]);
        fprintf(fp, "%10.8e,", magneticforce[i]);
        fprintf(fp, "%10.8e,", current[i]);
        fprintf(fp, "%10.8e,", flux[i]);
        fprintf(fp, "%d,", dynamicstepsarray[i]);
        fprintf(fp, "%10.8e,", time[i]);
        fprintf(fp, "; \n");
    }
    fprintf(fp, "];\n");

    fprintf(fp, "subplot(2,3,1);hold on;\n");
    fprintf(fp, "plot(results(:,1),results(:,2),'*-');\n");
    //fprintf(fp, "plot(results(:,1),results(:,3),'*-');\n");
    fprintf(fp, "title(\"%s\");\n\n", "displacement");

    fprintf(fp, "subplot(2,3,2);hold on;\n");
    fprintf(fp, "plot(results(:,1),results(:,3),'*-');\n");
    fprintf(fp, "title(\"%s\");\n\n", "velocities");

    fprintf(fp, "subplot(2,3,3);hold on;\n");
    fprintf(fp, "plot(results(:,1),results(:,4),'*-');\n");
    fprintf(fp, "title(\"%s\");\n\n", "accelerations");

    fprintf(fp, "subplot(2,3,4);hold on;\n");
    fprintf(fp, "plot(results(:,1),results(:,5),'*-');\n");
    fprintf(fp, "title(\"%s\");\n\n", "magforces");

    fprintf(fp, "subplot(2,3,5);hold on;\n");
    fprintf(fp, "plot(results(:,1),results(:,6),'*-');\n");
    //fprintf(fp, "plot(results(:,1),results(:,9),'*-');\n");
    fprintf(fp, "title(\"%s\");\n\n", "ICoil");

    fprintf(fp, "subplot(2,3,6);hold on;\n");
    fprintf(fp, "plot(results(:,1),results(:,7),'*-');\n");
    fprintf(fp, "title(\"%s\");\n", "PhiCoil");

    fclose(fp);
    cudaFree(d_Flux);
    cudaFree(d_FluxArray);
}

void solveWeakCouple(int step)
{
    double i_tmp, flux_tmp, L_tmp, dfluxdt, f, dfdi;
    if (step == 1) {
        i_tmp = 0.01;
    }
    else {
        i_tmp = current[step - 1];
    }

    int maxstep = 20;

    double* a, * b, error, * d_error;
    cudaMalloc(&a, m_num_nodes * sizeof(double));
    cudaMalloc(&b, m_num_nodes * sizeof(double));
    cudaMalloc(&d_error, sizeof(double));
    for (int i = 0; i < maxstep; ++i) {
        error = 1;
        //updateLoad << <CudaBlckNum, CudaThrdNum >> > (m_num_triele, d_mp_triele, i_tmp);
        for (int n = 0; n < m_num_triele; ++n) {
            double Jor = i_tmp * 2400 / 0.0147 / 0.011687;
            if (mp_triele[n].domain == 7) {
                mp_triele[n].J = Jor;
            }
        }
        cudaMemcpy(d_mp_triele, mp_triele, m_num_triele * sizeof(CTriElement), cudaMemcpyHostToDevice);

        cudaDeviceSynchronize();
        for (int iter = 0; iter < 30000; ++iter) {
            UpdateSolutiontoA1 << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_triele, d_mp_node);
            //calculateGlobalError << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_triele, d_mp_node, a, b, d_error);
            //cudaMemcpy(&error, d_error, sizeof(double), cudaMemcpyDeviceToHost);
            //cudaDeviceSynchronize();
            if (iter % 100 == 0) {
                calculateGlobalError << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_triele, d_mp_node, a, b, d_error);
                cudaMemcpy(&error, d_error, sizeof(double), cudaMemcpyDeviceToHost);
                printf("iter: %d, error: %.20f\n", iter, error);
            }

            //printf("iter: %d, error: %.20f\n", iter, error);
            UpdateAttoAtold << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_node);
            if (error < 1e-5) {
                staticsteps = iter;
                break;
            }
        }
        cudaDeviceSynchronize();
        dynamicsteps += staticsteps;
        //solveFlux << <CudaBlckNum, CudaThrdNum >> > (m_num_triele, d_mp_triele, d_mp_node, d_FluxArray, d_Flux);
        //cudaDeviceSynchronize();
        //cudaMemcpy(&flux_tmp, &d_FluxArray[0], sizeof(double), cudaMemcpyDeviceToHost);
        //flux_tmp = *d_Flux;

        cudaMemcpy(mp_node, d_mp_node, m_num_nodes * sizeof(CNode), cudaMemcpyDeviceToHost);
        //for (int n = 0; n < m_num_nodes; ++n) {
        //    cout << "n: " << n << "At: " << mp_node[n].At << endl;
        //}
        /*solve flux*/
        flux_tmp = 0;
        CTriElement triele;
        int n0, n1, n2;
        for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
            triele = mp_triele[i_tri];
            if (triele.domain != 7)
                continue;
            n0 = triele.n[0];
            n1 = triele.n[1];
            n2 = triele.n[2];
            //cout << triele.material->getFEMCoil().tau;
            flux_tmp += (2 * 3.1415927 * 13969821.69 * triele.area * (mp_node[n0].At + mp_node[n1].At + mp_node[n2].At) / 3);
        }
        /*solve flux finish*/
        L_tmp = flux_tmp / i_tmp;
        dfluxdt = (flux_tmp - flux[step - 1]) / h;
        f = i_tmp * R + dfluxdt - U;
        dfdi = R + L_tmp / h;
        i_tmp = i_tmp - f / dfdi;
        cout << "flux_tmp: " << flux_tmp << ", i_tmp: " << i_tmp << ", error: " << abs(f / dfdi / i_tmp) << endl;
        if (abs(f / dfdi / i_tmp) < 1e-6)
            break;
    }
    current[step] = i_tmp;
    flux[step] = flux_tmp;
    cudaFree(a);
    cudaFree(b);
    cudaFree(d_error);
}

double solveSpringForce(double pos)
{
    if (pos > position[4 - 1]) {
        int  i = 4 - 1;
        double k = (force[i] - force[i - 1]) / (position[i] - position[i - 1]);
        return force[i] + k * (pos - position[i]);
    }
    else if (pos < position[0]) {
        return force[0];
    }
    else
    {
        for (int i = 0; i < 4 - 1; i++) {
            if (pos >= position[i] && pos <= position[i + 1]) {
                double k = (force[i + 1] - force[i]) / (position[i + 1] - position[i]);
                return force[i] + k * ((pos - position[i]));
            }
        }
    }
}

double solveMagneticForce()
{
    vector<bool> deformednode(m_num_nodes, false);
    vector<bool> movenode(m_num_nodes, false);
    vector<int> deformelement;
    vector<int> moveelement;

    //遍历全部单元，标记节点
    for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
        if (mp_triele[i_tri].domain == 3) {
            deformelement.push_back(i_tri);
            //把三个节点标记为形变区域节点
            for (int i = 0; i < 3; ++i) {
                deformednode[mp_triele[i_tri].n[i]] = true;
            }
        }
    }

    for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
        if (mp_triele[i_tri].domain == 4) {
            moveelement.push_back(i_tri);
            //把三个节点标记为运动区域节点
            for (int i = 0; i < 3; ++i) {
                movenode[mp_triele[i_tri].n[i]] = true;
            }
        }
    }

    //将全部节点的电磁力初始化为0
    for (int n = 0; n < m_num_nodes; ++n) {
        mp_node[n].NodeForcex = 0;
        mp_node[n].NodeForcey = 0;
    }

    //遍历全部形变区域单元，如果单元内节点也处在运动单元上,
    //计算相应的单元电磁力和节点电磁力
    //师兄代码里的A是At吗？？？
    for (auto i_tri : deformelement) {
        CTriElement triele = mp_triele[i_tri];
        double rc = triele.rc;
        double area = triele.area;
        double mu = 4 * 3.1415927 * 1e-7;
        double QA = triele.Q[0] * mp_node[triele.n[0]].At + triele.Q[1] * mp_node[triele.n[1]].At + triele.Q[2] * mp_node[triele.n[2]].At;
        double RA = triele.R[0] * mp_node[triele.n[0]].At + triele.R[1] * mp_node[triele.n[1]].At + triele.R[2] * mp_node[triele.n[2]].At;
        double tmp = 3.1415927 / 4 / mu;	//为什么要除以4???
        double dvdB2 = 0;
        if (dvdB2 != 0) {
            cout << "i_tri: " << i_tri << endl;
        }
        double tmp1 = 3.1415927 / 4 * (QA * QA + RA * RA) / area / rc * dvdB2 / 4;
        double dQdy[3][3];
        double dRdx[3][3];
        double ac = area * rc;
        double ac2 = ac * ac;
        double ac3 = ac2 * ac;

        dQdy[0][0] = 0; dQdy[0][1] = -1; dQdy[0][2] = 1;
        dQdy[1][0] = 1; dQdy[1][1] = 0; dQdy[1][2] = -1;
        dQdy[2][0] = -1; dQdy[2][1] = 1; dQdy[2][2] = 0;

        dRdx[0][0] = 0; dRdx[0][1] = 1; dRdx[0][2] = -1;
        dRdx[1][0] = -1; dRdx[1][1] = 0; dRdx[1][2] = 1;
        dRdx[2][0] = 1; dRdx[2][1] = -1; dRdx[2][2] = 0;

        double dSdx, dSdy;
        double dQAdy, dRAdx;
        double dBxdx, dBydx, dBxdy, dBydy;
        for (int i = 0; i < 3; ++i) {
            int n = mp_triele[i_tri].n[i];
            if (movenode[n] == true) {
                //x方向
                dSdx = 0.5 * (dRdx[i][2] * triele.Q[1] - dRdx[i][1] * triele.Q[2]);
                dRAdx = (dRdx[i][0] * mp_node[triele.n[0]].At + dRdx[i][1] * mp_node[triele.n[1]].At + dRdx[i][2] * mp_node[triele.n[2]].At);
                mp_node[n].NodeForcex -= tmp * (2 * RA / ac * dRAdx);
                mp_node[n].NodeForcex -= tmp * ((QA * QA + RA * RA) * (-1 / ac2 * (area / 3 + rc * dSdx)));
                mp_node[n].NodeForcex -= tmp1 * (2 * RA / ac2 * dRAdx);
                mp_node[n].NodeForcex -= tmp1 * ((QA * QA + RA * RA) * (-2 / ac3 * (area / 3 + rc * dSdx)));
                //y方向
                dSdy = 0.5 * (dQdy[i][1] * triele.R[2] - dQdy[i][2] * triele.R[1]);
                dQAdy = (dQdy[i][0] * mp_node[triele.n[0]].At + dQdy[i][1] * mp_node[triele.n[1]].At + dQdy[i][2] * mp_node[triele.n[2]].At);
                mp_node[n].NodeForcey -= tmp * ((QA * QA + RA * RA) * (-1) / ac2 * rc * dSdy);	//为什么会有负号？
                mp_node[n].NodeForcey -= tmp * (2 * QA / ac * dQAdy);
                mp_node[n].NodeForcey -= tmp1 * ((QA * QA + RA * RA) * (-2) / ac3 * rc * dSdy);
                mp_node[n].NodeForcey -= tmp1 * (2 * QA / ac2 * dQAdy);
            }
        }
    }

    double Fy = 0;
    for (int n = 0; n < m_num_nodes; ++n) {
        if (movenode[n] == true) {
            Fy += mp_node[n].NodeForcey;
        }
    }
    return Fy;
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
    int NRCount = 10;
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


                if (d_MyElem[ID].linearflag == true)	
                {
                    Vt = 1.0 / (4 * 3.14159 * 1e-7) / d_MyElem[ID].xdot;
                    J0 = Vt * d_MyElem[ID].C[nodenumber][nodenumber];
                    k0 = Js - Vt * RHSContri;
                    //V = 1.0 / d_MyElem[ID].material->getMu(0);
                    //J0 = V * d_MyElem[ID].C[nodenumber][nodenumber];
                    //k0 = Js - V * RHSContri;
                }

                else	
                {
                    B2 = -1 / d_MyElem[ID].area * (d_MyElem[ID].C[0][1] * (AtLocal[0] - AtLocal[1]) * (AtLocal[0] - AtLocal[1]) + d_MyElem[ID].C[1][2] * (AtLocal[1] - AtLocal[2]) * (AtLocal[1] - AtLocal[2]) + d_MyElem[ID].C[0][2] * (AtLocal[0] - AtLocal[2]) * (AtLocal[0] - AtLocal[2]));
                    B = sqrt(B2) / d_MyElem[ID].xdot;
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

                    B2A = -1 / d_MyElem[ID].area * (2 * d_MyElem[ID].C[nodenumber][0] * (At - AtLocal[0]) + 2 * d_MyElem[ID].C[nodenumber][1] * (At - AtLocal[1]) + 2 * d_MyElem[ID].C[nodenumber][2] * (At - AtLocal[2])) / d_MyElem[ID].xdot;
                    J0 = B2A * VB2 * RHSContri / d_MyElem[ID].xdot / d_MyElem[ID].xdot + Vt * d_MyElem[ID].C[nodenumber][nodenumber];
                    k0 = Js - Vt * RHSContri;

                    //B2 = -1 / d_MyElem[ID].area * (d_MyElem[ID].C[0][1] * (AtLocal[0] - AtLocal[1]) * (AtLocal[0] - AtLocal[1]) + d_MyElem[ID].C[1][2] * (AtLocal[1] - AtLocal[2]) * (AtLocal[1] - AtLocal[2]) + d_MyElem[ID].C[0][2] * (AtLocal[0] - AtLocal[2]) * (AtLocal[0] - AtLocal[2]));
                    //B = sqrt(B2);	//����B�����ڴ���������
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
            if (NRerror > 1e-6) {
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

    ////�ж�ȫ��������
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
    //    //��һ��Ӱ��Ч��
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
    a[n] = 0;
    b[n] = 0;
    __syncthreads();

    a[n] = (d_MyNode[n].At - d_MyNode[n].At_old) * (d_MyNode[n].At - d_MyNode[n].At_old);
    b[n] = d_MyNode[n].At * d_MyNode[n].At;
    //d_error = 0;
    __syncthreads();

    //a��b��Լ���
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
        *d_error = sqrt(a[0]) / sqrt(b[0]);
        //printf("a: %.30f, b: %.30f\n", a[0], b[0]);
    }
}

__global__ void updateLoad(int m_num_triele, CTriElement* d_MyElem, double current)
{
    int n = blockDim.x * blockIdx.x + threadIdx.x;
    if (n >= m_num_triele) {
        return;
    }

    double Jor = current * 2400 / 0.0147 / 0.011687;
    if (d_MyElem[n].domain == 7) {
        d_MyElem[n].J = Jor;
    }
}

__global__ void solveFlux(int m_num_triele, CTriElement* d_MyElem, CNode* d_MyNode, double* d_FluxArray, double* d_Flux)
{
    int i_tri = blockDim.x * blockIdx.x + threadIdx.x;
    if (i_tri >= m_num_triele) {
        return;
    }

    if (d_MyElem[i_tri].domain == 7) {
        int n0 = d_MyElem[i_tri].n[0];
        int n1 = d_MyElem[i_tri].n[1];
        int n2 = d_MyElem[i_tri].n[2];
        d_FluxArray[i_tri] = (d_MyNode[n0].At + d_MyNode[n1].At + d_MyNode[n2].At) / 3;
        d_FluxArray[i_tri] *= (2 * 3.1415927 * 13969821.69 * d_MyElem[i_tri].area);

    }
    else {
        d_FluxArray[i_tri] = 0;
    }
    __syncthreads();

    int leng = m_num_triele;
    for (int i = m_num_triele / 2.0 + 0.5; i > 1; i = i / 2.0 + 0.5) {
        if (i_tri < i)
        {

            if (i_tri + i < leng)
            {
                d_FluxArray[i_tri] += d_FluxArray[i_tri + i];
            }
        }
        __syncthreads();
        leng = leng / 2.0 + 0.5;
    }
    __syncthreads();
    if (i_tri == 0) {
        d_FluxArray[0] = d_FluxArray[0] + d_FluxArray[1];
        *d_Flux = d_FluxArray[0];
        printf("d_FluxArray: %f\n", d_FluxArray[0]);
    }
}

__device__ double getV(double B) {
    //if (B <= 0.6) {
    //    return 500;
    //}
    //else {
    //    return 500 + 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
    //}

    double slope, H, b;
    if (B < 1e-3)  return Hdata[1] / Bdata[1];
    getkHb(B, &slope, &H, &b);
    if (B / H < 3.1415927 * 4e-7)  return 1 / (3.1415927 * 4e-7);
    return H / B;

    //if (B <= 0.6) {
    //    double H = 663.146 * B;
    //    return 663.146;
    //}
    //else if (B > 0.6) {
    //    double H = 663.146 * B + 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6);
    //    return H / B;
    //    ////return 0.002;
    //}
}

__device__ double getdVdB2(double B) {
    //if (B <= 0.6) {
    //    return 0;
    //}
    //else {
    //    return (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
    //}

    double slope, H, b;
    if (B < 1e-9) return 0;
    getkHb(B, &slope, &H, &b);
    return -b / (B * B * B * 2);

    //if (B <= 0.6) {
    //    return 0;
    //}
    //else if (B > 0.6) {
    //    double dvdB2 = (18000 * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6) * B - 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
    //    return dvdB2;
    //}
}

__device__ double getkHb(double B, double* k, double* H, double* b) {
    if (B >= Bdata[BHpoints - 1]) {
        int  i = BHpoints - 1;
        (*k) = (Hdata[i] - Hdata[i - 1]) / (Bdata[i] - Bdata[i - 1]);
        (*b) = Hdata[i - 1] - (*k) * Bdata[i - 1];
    }
    else if (B < Bdata[0]) {
        (*k) = (Hdata[1] - Hdata[0]) / (Bdata[1] - Bdata[0]);
        (*b) = Hdata[0] - (*k) * Bdata[0];
    }
    else {
        for (int i = 0; i < BHpoints - 1; i++) {
            if (B >= Bdata[i] && B <= Bdata[i + 1]) {
                (*k) = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
                (*b) = Hdata[i] - (*k) * Bdata[i];
                break;
            }
        }
    }
    (*H) = (*k) * B + (*b);
}