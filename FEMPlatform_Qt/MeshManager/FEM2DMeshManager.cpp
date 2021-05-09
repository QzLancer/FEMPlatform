#include "FEM2DMeshManager.h"


FEM2DMeshManager::~FEM2DMeshManager()
{

}

void FEM2DMeshManager::readMeshFile(string meshfile)
{
	//获取分网文件拓展名
    if (meshfile.empty()) {
        meshfile = this->meshfile;
    }
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
    //分网文件为.msh的情况
    else if (extension == ".msh") {
        read2DMsh(meshfile);
    }
	else {
		std::cout << "Error: Can't read this meshfile type!" << endl;
        exit(0);
	}


}

void FEM2DMeshManager::read2DMphtxt(string meshfile)
{
    if (mp_node != nullptr) delete[] mp_node;
    if (mp_vtxele != nullptr) delete[] mp_vtxele;
    if (mp_edgele != nullptr) delete[] mp_edgele;
    if (mp_triele != nullptr) delete[] mp_triele;
    cout << meshfile << endl;
    FILE* fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, meshfile.c_str(), "r");
    if (!fp) {
        std::cout << "Error: opening file!" << endl;
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

void FEM2DMeshManager::read2DMsh(string meshfile)
{
    if (mp_node != nullptr) delete[] mp_node;
    if (mp_vtxele != nullptr) delete[] mp_vtxele;
    if (mp_edgele != nullptr) delete[] mp_edgele;
    if (mp_triele != nullptr) delete[] mp_triele;
    cout << "Read Gmsh mesh file: " << meshfile << endl;
    FILE* fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, meshfile.c_str(), "r");
    if (!fp) {
        std::cout << "Error: opening file!" << endl;
        exit(0);
    }

    long offset; //文本读取位置标记
    int ele_number;
    int ele_type;
    int number_of_tags;
    int physic_tag;
    int geometry_tag;
    char* chtmp;
    m_num_vtxele = 0;
    m_num_edgele = 0;
    m_num_triele = 0;

    while (!feof(fp)) {
        fgets(ch, 256, fp);

        if (strstr(ch, "$MeshFormat")) {
            double version;
            int file_type;
            int data_size;
            if (fscanf_s(fp, "%lf %d %d\n", &version, &file_type, &data_size) != 3) {
                printf("error reading format");
                return;
            }
            else {
                if (version > 2.2) {
                    printf("Can only open gmsh version 2.2 format");
                    return;
                }
            }
            fgets(ch, 256, fp);
            if (!strstr(ch, "$EndMeshFormat")) {
                printf("$MeshFormat section should end to string $EndMeshFormat:\n%s\n", ch);
            }
        }
        else if (strstr(ch, "$Nodes")) {
            if (fscanf_s(fp, "%d\n", &m_num_nodes) != 1) {
                return;
            }
            else {
                /** 读取节点坐标 **/
                mp_node = new CNode[m_num_nodes];
                int index;
                for (int i = 0; i < m_num_nodes; ++i) {
                    fscanf_s(fp, "%d %lf %lf %lf\n", &index, &mp_node[i].x, &mp_node[i].y, &mp_node[i].z);
                    //qDebug()<<index<<pmeshnode[i].x<<pmeshnode[i].y<<pmeshnode[i].z;
                }
            }
            fgets(ch, 256, fp);
            if (!strstr(ch, "$EndNodes")) {
                printf("$Node section should end to string $EndNodes:\n%s\n", ch);
            }
        }
        else if (strstr(ch, "$Elements")) {
            if (fscanf_s(fp, "%d\n", &m_num_ele) != 1) {
                return;
            }
            else {
                offset = ftell(fp);
                for (int i = 0; i < m_num_ele; ++i) {
                    chtmp = fgets(ch, 256, fp);
                    ele_number = next_int(&chtmp);
                    ele_type = next_int(&chtmp);
                    number_of_tags = next_int(&chtmp);
                    physic_tag = next_int(&chtmp);
                    geometry_tag = next_int(&chtmp);
                    switch (ele_type) {
                    case 15:/** point **/
                        m_num_vtxele++;
                        break;
                    case 1:/** line **/
                        m_num_edgele++;
                        break;
                    case 2:/** triangle **/
                        m_num_triele++;
                        break;
                    }
                }

                std::cout << m_num_nodes << " number of nodes." << endl;
                std::cout << m_num_vtxele << " number of vertex elements." << endl;
                std::cout << m_num_edgele << " number of edge elements." << endl;
                std::cout << m_num_triele << " number of triangle elements." << endl;

                mp_vtxele = new CVtxElement[m_num_vtxele];
                mp_edgele = new CEdgElement[m_num_edgele];
                mp_triele = new CTriElement[m_num_triele];

                //回到offset处，将信息写入数组
                fseek(fp, offset, SEEK_SET);
                int vtxcnt = 0, edgcnt = 0, tricnt = 0;
                for (int i = 0; i < m_num_ele; ++i) {
                    chtmp = fgets(ch, 256, fp);
                    ele_number = next_int(&chtmp);
                    ele_type = next_int(&chtmp);
                    number_of_tags = next_int(&chtmp);
                    physic_tag = next_int(&chtmp);
                    geometry_tag = next_int(&chtmp);
                    switch (ele_type)
                    {
                    case 15:
                        mp_vtxele[vtxcnt].domain = geometry_tag;
                        mp_vtxele[vtxcnt].n = next_int(&chtmp) - 1;
                        vtxcnt++;
                        break;
                    case 1:
                        mp_edgele[edgcnt].domain = geometry_tag;
                        for (int j = 0; j < 2; ++j) {
                            mp_edgele[edgcnt].n[j] = next_int(&chtmp) - 1;
                        }
                        edgcnt++;
                        break;
                    case 2:
                        mp_triele[tricnt].domain = geometry_tag;
                        for (int j = 0; j < 3; ++j) {
                            mp_triele[tricnt].n[j] = next_int(&chtmp) - 1;
                        }
                        tricnt++;
                        break;
                    default:
                        break;
                    }
                }

            }
            fgets(ch, 256, fp);
            if (!strstr(ch, "$EndElements")) {
                printf("$Element section should end to string $EndElements:\n%s\n", ch);
            }
        }
    }




    fclose(fp);


}
