#include "FEM2DNDDRCUDASolver.cuh"

#include <time.h>
#include <stdio.h>

void FEM2DNDDRCUDASolver::solve()
{
	GPUInitialMallocCopy();
	//这些初始化部分能不能改成并行呢？
	makeTrangle();	
	processBoundaryCondition();	//似乎不要用到边界点排序这一过程
	//clock_t start, end;
	//start = clock();
	processMaterial();	//执行速度非常慢
	//end = clock();
	//cout << "process material time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
	processLoad();
	//处理NDDR的节点信息
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
	for (int iter = 0; iter < maxitersteps; ++iter) {
		nodeAnalysis <<<CudaBlckNum, CudaThrdNum >>> (m_num_nodes, mp_node, mp_triele);
		cudaDeviceSynchronize();
		//判断全局收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
			b += mp_node[i].At * mp_node[i].At;
		}
		error = sqrt(a) / sqrt(b);
	
		if ((iter + 1) % 100 == 0) {
			cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		}
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].At_old = mp_node[i].At;
			}
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			return;
		}

		
	}
}

FEM2DNDDRCUDASolver::~FEM2DNDDRCUDASolver()
{
	if (mp_triele != nullptr) {
		cudaFree(mp_triele);
		mp_triele = nullptr;
	}
	if (mp_edgele != nullptr) {
		cudaFree(mp_edgele);
		mp_edgele = nullptr;
	}
	if (mp_vtxele != nullptr) {
		cudaFree(mp_vtxele);
		mp_vtxele = nullptr;
	}
	if (mp_node != nullptr) {
		cudaFree(mp_node);
		mp_node = nullptr;
	}
}

void FEM2DNDDRCUDASolver::setNodes(const int _numofnodes, CNode* const _nodes)
{
	cout << "FEM2DNDDRCUDASolver::setNodes\n";
	m_num_nodes = _numofnodes;
	cudaMallocManaged((void**)&mp_node, m_num_nodes * sizeof(CNode));
	memcpy(mp_node, _nodes, m_num_nodes * sizeof(CNode));
}

void FEM2DNDDRCUDASolver::setVtxElements(const int _numofvtx, CVtxElement* const _vtxele)
{
	m_num_vtxele = _numofvtx;
	cudaMallocManaged((void**)&mp_vtxele, m_num_vtxele * sizeof(CVtxElement));
	memcpy(mp_vtxele, _vtxele, m_num_vtxele * sizeof(CVtxElement));
}

void FEM2DNDDRCUDASolver::setEdgElements(const int _numofedg, CEdgElement* const _edgele)
{
	m_num_edgele = _numofedg;
	cudaMallocManaged((void**)&mp_edgele, m_num_edgele * sizeof(CEdgElement));
	memcpy(mp_edgele, _edgele, m_num_edgele * sizeof(CEdgElement));
}

void FEM2DNDDRCUDASolver::setTriElements(const int _numoftri, CTriElement* const _triele)
{
	cout << "FEM2DNDDRCUDASolver::setTriElements(const int _numoftri, CTriElement* const _triele)\n";
	m_num_triele = _numoftri;
	cudaMallocManaged((void**)&mp_triele, m_num_triele * sizeof(CTriElement));
	memcpy(mp_triele, _triele, m_num_triele * sizeof(CTriElement));
}

void FEM2DNDDRCUDASolver::processMaterial()
{
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int domain = mp_triele[i_tri].domain;
		FEMMaterial* material = materialmap[domain];
		cudaMallocManaged(&mp_triele[i_tri].material, sizeof(FEMMaterial));
		mp_triele[i_tri].material->GPUCopy(*material);
	}
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

}

__global__ void nodeAnalysis(int d_m_num_nodes, CNode* d_mp_node, CTriElement* d_mp_triele)
{
	int n = threadIdx.x + blockIdx.x * blockDim.x;
	if (n >= d_m_num_nodes)
		return;
	if (d_mp_node[n].bdr == 1)
		return;
	//节点内部迭代过程
	int maxNRitersteps = 100;
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
							sigmai += triele.C[i][m] * d_mp_node[triele.n[i]].At_old;
						}
					}
					for (int j = 0; j < 3; ++j) {
						for (int m = 0; m < 3; ++m) {
							if (m == nodenumber) {
								sigmaj += triele.C[j][m] * Ati;
							}
							else {
								sigmaj += triele.C[j][m] * d_mp_node[triele.n[i]].At_old;
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
