#include "FEM2DNDDRCUDASolver.cuh"

#include <time.h>
#include <stdio.h>

__managed__ double error = 0;

void FEM2DNDDRCUDASolver::solve()
{
	//CPU处理初始条件，然后数据拷贝到GPU中
	processBoundaryCondition();
	processLoad();
	processNDDRNode();
	makeTrangle();
	processMaterial();

	GPUInitialMallocCopy();

	//makeTrangleinDevice << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_node, d_mp_triele);
	//cudaDeviceSynchronize();


	double* a, * b;
	cudaMalloc(&a, m_num_nodes * sizeof(double));
	cudaMalloc(&b, m_num_nodes * sizeof(double));
	error = 0;
	for (int iter = 0; iter < maxitersteps; ++iter) {
		nodeAnalysis << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_node, d_mp_triele);
		calculateGlobalError << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_node, d_mp_triele, a, b);
		cudaDeviceSynchronize();

		if ((iter + 1) % 100 == 0) {
			cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		}
		if (error > maxerror) {
			copyAttoAtold << <CudaBlckNum, CudaThrdNum >> > (m_num_nodes, d_mp_node);
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			break;
		}
	}

	cudaMemcpy(mp_node, d_mp_node, m_num_nodes * sizeof(CNode), cudaMemcpyDeviceToHost);

	cudaFree(b);
	cudaFree(a);
	GPUFree();
}

FEM2DNDDRCUDASolver::~FEM2DNDDRCUDASolver()
{
	//if (mp_triele != nullptr) {
	//	cudaFree(mp_triele);
	//	mp_triele = nullptr;
	//}
	//if (mp_edgele != nullptr) {
	//	cudaFree(mp_edgele);
	//	mp_edgele = nullptr;
	//}
	//if (mp_vtxele != nullptr) {
	//	cudaFree(mp_vtxele);
	//	mp_vtxele = nullptr;
	//}
	//if (mp_node != nullptr) {
	//	cudaFree(mp_node);
	//	mp_node = nullptr;
	//}
}

//void FEM2DNDDRCUDASolver::setNodes(const int _numofnodes, CNode* const _nodes)
//{
//	cout << "FEM2DNDDRCUDASolver::setNodes\n";
//	m_num_nodes = _numofnodes;
//	//cudaMallocManaged((void**)&mp_node, m_num_nodes * sizeof(CNode));
//	//memcpy(mp_node, _nodes, m_num_nodes * sizeof(CNode));
//	cudaMalloc(&mp_node, m_num_nodes * sizeof(CNode));
//	cudaMemcpy(mp_node, _nodes, m_num_nodes * sizeof(CNode), cudaMemcpyHostToDevice);
//}
//
//void FEM2DNDDRCUDASolver::setVtxElements(const int _numofvtx, CVtxElement* const _vtxele)
//{
//	m_num_vtxele = _numofvtx;
//	//cudaMallocManaged((void**)&mp_vtxele, m_num_vtxele * sizeof(CVtxElement));
//	//memcpy(mp_vtxele, _vtxele, m_num_vtxele * sizeof(CVtxElement));
	//cudaMalloc(&mp_vtxele, m_num_vtxele * sizeof(CVtxElement));
	//cudaMemcpy(mp_vtxele, _vtxele, m_num_vtxele * sizeof(CVtxElement), cudaMemcpyHostToDevice);
//}
//
//void FEM2DNDDRCUDASolver::setEdgElements(const int _numofedg, CEdgElement* const _edgele)
//{
//	m_num_edgele = _numofedg;
//	//cudaMallocManaged((void**)&mp_edgele, m_num_edgele * sizeof(CEdgElement));
//	//memcpy(mp_edgele, _edgele, m_num_edgele * sizeof(CEdgElement));
//	cudaMalloc((void**)&mp_edgele, m_num_edgele * sizeof(CEdgElement));
//	cudaMemcpy(mp_edgele, _edgele, m_num_edgele * sizeof(CEdgElement), cudaMemcpyHostToDevice);
//}
//
//void FEM2DNDDRCUDASolver::setTriElements(const int _numoftri, CTriElement* const _triele)
//{
//	cout << "FEM2DNDDRCUDASolver::setTriElements(const int _numoftri, CTriElement* const _triele)\n";
//	m_num_triele = _numoftri;
//	//cudaMallocManaged((void**)&mp_triele, m_num_triele * sizeof(CTriElement));
//	//memcpy(mp_triele, _triele, m_num_triele * sizeof(CTriElement));
	//cudaMalloc((void**)&mp_triele, m_num_triele * sizeof(CTriElement));
	//cudaMemcpy(mp_triele, _triele, m_num_triele * sizeof(CTriElement), cudaMemcpyHostToDevice);
//}

//在GPU中创建BH数组
//在GPU中创建材料数组，材料的BH曲线指向GPU
//目前是根据domain的大小来分配材料空间，还有进一步优化的可能
void FEM2DNDDRCUDASolver::processMaterial()
{
	//for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
	//	int domain = mp_triele[i_tri].domain;
	//	cout << "domain: " << domain << endl;
	//	FEMMaterial* material = materialmap[domain];
	//	cudaMallocManaged(&mp_triele[i_tri].material, sizeof(FEMMaterial));
	//	mp_triele[i_tri].material->GPUCopy(*material);
	//}
	int domainsize = 0;
	for (auto mat = materialmap.begin(); mat != materialmap.end(); mat++ ) {
		if (domainsize < mat->first) {
			domainsize = mat->first;
		}
	}

	//在CPU中创建材料数组
	cout << "domainsize: " << domainsize << endl;
	materialarray = new FEMMaterial[domainsize];

	for (int i = 0; i < domainsize; ++i) {
		materialarray[i].GPUCopy(*materialmap[i + 1]);
	}

	//将材料数组拷贝到GPU中
	cudaMalloc(&d_materialarray, domainsize * sizeof(FEMMaterial));
	cudaMemcpy(d_materialarray, materialarray, domainsize * sizeof(FEMMaterial), cudaMemcpyHostToDevice);

	//delete[] materialarray;
	
	//将material绑定到TriElement中
	//此为GPU版本，考虑到processLoad和Boundary的通用性，三者都放在GPU内存拷贝之前,因此先使用CPU版本。
	//GPU版本中，GPUInitialMallocCopy()放在processMaterial()之前
	//GPU版本
	//assignMattoTriEle << <CudaBlckNum, CudaThrdNum >> > (m_num_triele, d_mp_triele, d_materialarray);
	//cudaDeviceSynchronize();
 
	//此为CPU版本
	//CPU版本中，GPUInitialMallocCopy()放在processMaterial()之后
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int matid = mp_triele[i_tri].domain - 1;
		mp_triele[i_tri].material = &d_materialarray[matid];
	}

	//materialArrayTest << <CudaBlckNum, CudaThrdNum >> > (domainsize, d_materialarray);
	//cudaDeviceSynchronize();
}

void FEM2DNDDRCUDASolver::GPUFree()
{
	cudaFree(d_mp_node);
	cudaFree(d_mp_triele);
	cudaFree(d_materialarray);
	delete[] materialarray;
}

void FEM2DNDDRCUDASolver::GPUInitialMallocCopy()
{
	//CUDA initialize
	int num_devices, device;
	cudaGetDeviceCount(&num_devices);
	std::cout << "Number of device: " << num_devices << endl;
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
	for (int i = 0; i < num_devices; i++) {
		cudaDeviceProp prop;
		if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
			if (prop.major >= 1) {
				break;
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

void FEM2DNDDRCUDASolver::processNDDRNode()
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

__global__ void nodeAnalysis(int d_m_num_nodes, CNode* d_mp_node, CTriElement* d_mp_triele)
{
	int n = threadIdx.x + blockIdx.x * blockDim.x;
	if (n >= d_m_num_nodes)
		return;
	if (d_mp_node[n].bdr == 1)
		return;
	//节点内部迭代过程
	int maxNRitersteps = 1;
	double Ati = 0;
	for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
		double S = 0, F = 0;
		double J = 0, Fj = 0;
		//装配过程
		for (int k = 0; k < d_mp_node[n].NumberofNeighbourElement; ++k) {
			int i_tri = d_mp_node[n].NeighbourElementId[k];
			CTriElement triele = d_mp_triele[i_tri];
			int nodenumber = d_mp_node[n].NeighbourElementNumber[k];
			double mu = triele.material->getMuinDevice(d_mp_triele[i_tri].B);
			double mut = triele.material->getMuinDevice(d_mp_triele[i_tri].B) * triele.xdot;
			//printf("nodeid: %d, nodenumber: %d, mut: %f\n",n, nodenumber, mut);
			//printf("triele.j: %f\n", triele.J);
			//处理线性单元
			if (triele.material->getLinearFlaginDevice() == true) {
				for (int i = 0; i < 3; ++i) {
					double Se = triele.C[nodenumber][i] / mut;
					if (nodenumber == i) {
						S += Se;
						F += triele.J * triele.area / 3;
					}
					else {
						F -= Se * d_mp_node[triele.n[i]].At_old;
					}
				}
			}
			//处理非线性单元
			else {
				double dvdb, dvdbt, Bt, sigmai = 0, sigmaj = 0;
				dvdb = triele.material->getdvdBinDevice(d_mp_triele[i_tri].B);
				dvdbt = dvdb / triele.xdot / triele.xdot;
				Bt = d_mp_triele[i_tri].B * triele.xdot;
				for (int i = 0; i < 3; ++i) {
					for (int m = 0; m < 3; ++m) {
						if (m == nodenumber) {
							sigmai += triele.C[i][m] * Ati;
						}
						else {
							sigmai += triele.C[i][m] * d_mp_node[triele.n[m]].At_old;
						}
					}
					for (int j = 0; j < 3; ++j) {
						for (int m = 0; m < 3; ++m) {
							if (m == nodenumber) {
								sigmaj += triele.C[j][m] * Ati;
							}
							else {
								sigmaj += triele.C[j][m] * d_mp_node[triele.n[m]].At_old;
							}
						}
					}
				}
				for (int i = 0; i < 3; ++i) {
					if (Bt != 0) {
						J = triele.C[nodenumber][i] / mut + sigmai * sigmaj / Bt / triele.area;
					}
					else {
						J = triele.C[nodenumber][i] / mut;
					}
					if (nodenumber == i) {
						S += J;
						F += (J - triele.C[nodenumber][i] / mut) * Ati;
					}
					else {
						F += (J - triele.C[nodenumber][i] / mut) * d_mp_node[triele.n[i]].At_old;
						F -= J * d_mp_node[triele.n[i]].At_old;
					}
				}
			}
		}
		//printf("NR_iter: %d\n", NRiter);
		//if (F != 0)
		//	printf("NR_iter: %d, nodeid: %d, S: %f, F: %f\n", NRiter, n, S, F);
		//Ati事实上不全部为0，但是无法输出更多小数点后位数
		Ati = F / S;
		//if (Ati != 0) {
		//	printf("NRiter: %d, nodeid: %d, S: %f, F: %f, Ati: %f\n", NRiter, n, S, F, Ati);
		//}
		//NR迭代收敛性判断
		double a = (Ati - d_mp_node[n].At) * (Ati - d_mp_node[n].At);
		double b = Ati * Ati;
		double NRerror = sqrtf(a) / sqrtf(b);
		//printf("Ati: %f, d_mp_node[n].At: %f, NRerror: %f\n", Ati, d_mp_node[n].At, NRerror);
		//__syncthreads();
		if (Ati == 0) {
			continue;
		}
		if (NRerror > 1e-5) {
			d_mp_node[n].At = Ati;
			d_mp_node[n].A = d_mp_node[n].At / d_mp_node[n].x;
			for (int i = 0; i < d_mp_node[n].NumberofNeighbourElement; ++i) {
				//updateB
				double bx = 0, by = 0;
				int i_tri = d_mp_node[n].NeighbourElementId[i];
				for (int j = 0; j < 3; ++j) {
					int n = d_mp_triele[i_tri].n[j];
					bx += d_mp_triele[i_tri].R[j] * d_mp_node[n].A;
					by += d_mp_triele[i_tri].Q[j] * d_mp_node[n].A;
				}
				bx = bx / 2 / d_mp_triele[i_tri].area;
				d_mp_triele[i_tri].Bx = bx;
				by = -by / 2 / d_mp_triele[i_tri].area;
				d_mp_triele[i_tri].By = by;
				d_mp_triele[i_tri].B = sqrtf(bx * bx + by * by);
			}
		}
		else {
			//printf("n:%d, NRiter: %d\n", n, NRiter);
			break;
		}
	}
}

__global__ void calculateGlobalError(int d_m_num_nodes, CNode* d_mp_node, CTriElement* d_mp_triele, double* a, double* b)
{
	int n = blockDim.x * blockIdx.x + threadIdx.x;
	if (n >= d_m_num_nodes) {
		return;
	}

	a[n] = (d_mp_node[n].At - d_mp_node[n].At_old) * (d_mp_node[n].At - d_mp_node[n].At_old);
	b[n] = d_mp_node[n].At * d_mp_node[n].At;
	__syncthreads();

	//a和b归约求和
	int leng = d_m_num_nodes;
	for (int i = d_m_num_nodes / 2.0 + 0.5; i > 1; i = i / 2.0 + 0.5) {
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
		error = sqrtf(a[0]) / sqrtf(b[0]);
	}
}

__global__ void copyAttoAtold(int d_m_num_nodes, CNode* d_mp_node)
{
	int n = blockDim.x * blockIdx.x + threadIdx.x;
	if (n >= d_m_num_nodes) {
		return;
	}

	d_mp_node[n].At_old = d_mp_node[n].At;
}

__global__ void makeTrangleinDevice(int numofTrangle, CNode* d_mp_node, CTriElement* d_mp_triele)
{
	//怎么判断保证thread数大于numofTrangle?
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index >= numofTrangle)
		return;

	int k, m, n;
	double q0, q1, q2, r0, r1, r2, area;
	k = d_mp_triele[index].n[0];
	m = d_mp_triele[index].n[1];
	n = d_mp_triele[index].n[2];

	q0 = d_mp_node[m].y - d_mp_node[n].y;
	d_mp_triele[index].Q[0] = q0;
	q1 = d_mp_node[n].y - d_mp_node[k].y;
	d_mp_triele[index].Q[1] = q1;
	q2 = d_mp_node[k].y - d_mp_node[m].y;
	d_mp_triele[index].Q[2] = q2;

	r0 = d_mp_node[n].x - d_mp_node[m].x;
	d_mp_triele[index].R[0] = r0;
	r1 = d_mp_node[k].x - d_mp_node[n].x;
	d_mp_triele[index].R[1] = r1;
	r2 = d_mp_node[m].x - d_mp_node[k].x;
	d_mp_triele[index].R[2] = r2;

	area = 0.5 * std::abs(q1 * r2 - r1 * q2);
	d_mp_triele[index].area = area;

	d_mp_triele[index].rc = (d_mp_node[k].x +
		d_mp_node[m].x +
		d_mp_node[n].x) / 3;
	d_mp_triele[index].zc = (d_mp_node[k].y +
		d_mp_node[m].y +
		d_mp_node[n].y) / 3;

	int flag = 0;
	for (int f = 0; f < 3; f++) {
		if (d_mp_node[d_mp_triele[index].n[f]].x < 1e-7) {
			flag++;
		}
	}

	//计算三角形重心半径
	if (flag == 2) {
		d_mp_triele[index].xdot = d_mp_triele[index].rc;
	}
	else {
		d_mp_triele[index].xdot = 1 / (d_mp_node[k].x + d_mp_node[m].x);
		d_mp_triele[index].xdot += 1 / (d_mp_node[k].x + d_mp_node[n].x);
		d_mp_triele[index].xdot += 1 / (d_mp_node[m].x + d_mp_node[n].x);
		d_mp_triele[index].xdot = 1.5 / d_mp_triele[index].xdot;
	}

	//计算一阶三角形轴对称单元系数矩阵
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			d_mp_triele[index].C[i][j] = ((d_mp_triele[index].R[i] * d_mp_triele[index].R[j] + d_mp_triele[index].Q[i] * d_mp_triele[index].Q[j])) / (4 * d_mp_triele[index].area);
		}
	}
}

__global__ void assignMattoTriEle(int numofTrangle, CTriElement* d_mp_triele, FEMMaterial* d_materialarray)
{
	int i_tri = blockDim.x * blockIdx.x + threadIdx.x;
	if (i_tri >= numofTrangle) {
		return;
	}

	int matid = d_mp_triele[i_tri].domain - 1;
	d_mp_triele[i_tri].material = &d_materialarray[matid];
}
