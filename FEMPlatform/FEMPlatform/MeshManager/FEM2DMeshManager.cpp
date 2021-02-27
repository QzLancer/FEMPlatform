#include "FEM2DMeshManager.h"


FEM2DGenerator::~FEM2DGenerator()
{

}

void FEM2DGenerator::readMeshFile(string meshfile)
{
	//获取分网文件拓展名
	string extension;
	for (int i = meshfile.size() - 1; i >= 0; --i) {
		if (meshfile[i] == '.') {
			extension = meshfile.substr(i, meshfile.size()-i);
			break;
		}
	}
	//分网文件为mphtxt的情况
	if (extension == ".mphtxt") {
		read2DMphtxt(meshfile);
	}
	else {
		std::cout << "Error: Can't read this meshfile type!" << endl;
        exit(0);
	}

}

void FEM2DGenerator::read2DMphtxt(string meshfile)
{
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
    mp_node = new C2DNode[m_num_nodes];
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
