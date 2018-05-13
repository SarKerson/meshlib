#pragma once
#include "MyMesh.h"
#include <eigen3/Eigen/Sparse>

using namespace MeshLib;


#define DISK 0
#define QUAD 1

class harmornicMap
{
private:
        
    double Energy;
    double e;
public:
/**
 * function of the boudary vertexs, [x0, y0, z0],[x1, y1, z1],...,[xn, yn, zn]
 */
    std::map<CMyVertex*, Eigen::Vector3d> g;    //all the boundary

    std::map<CMyVertex*, Eigen::Vector3d> f;    //all the result 

    std::map<CMyVertex*, Eigen::Vector3d> org;  //origin

public:

    harmornicMap(double e = 0.01): e(e), Energy(0) {}

    double k(CMyVertex * v, CMyVertex * w, CMyMesh & mesh)
    {
        CEdge * edge =  mesh.vertexEdge(v, w);
        Eigen::Vector3d v_v = org[v]; //v_v << (v->point())[0], (v->point())[1], (v->point())[2];
        Eigen::Vector3d v_w = org[w]; //v_w << (w->point())[0], (w->point())[1], (w->point())[2];
        assert(edge != NULL);
        if (v->boundary() && w->boundary()) {
            CHalfEdge * he = (edge->halfedge(0))->he_next();
            assert(he != NULL);
            CVertex * vk = he->vertex();
            Eigen::Vector3d v_k; v_k << (vk->point())[0], (vk->point())[1], (vk->point())[2];
            double mod1 = mod((v_v - v_k).cross(v_w - v_k));
            assert(mod1 != 0.0);
            double cot_k = (v_v - v_k).dot(v_w - v_k) / mod1;
            out << "mod: " << mod1 << "\n";
            
#if MESH_DEBUG
            std::cout << "cot: " << cot_k << "\n";
#endif
            return cot_k;
        }
        else {
            CHalfEdge * he0 = (edge->halfedge(0))->he_next(),
                      * he1 = (edge->halfedge(1))->he_next();
            assert(he0 != NULL && he1 != NULL);
            CVertex * vk = he0->vertex(),
                    * vl = he1->vertex();
            Eigen::Vector3d v_k = org[(CMyVertex*)vk]; //v_k << (vk->point())[0], (vk->point())[1], (vk->point())[2];
            Eigen::Vector3d v_l = org[(CMyVertex*)vl]; //v_l << (vl->point())[0], (vl->point())[1], (vl->point())[2];
            /* calculate the cot */
#if MESH_DEBUG
            std::cout << "vv - vk" << v_v - v_k << "\n"
                      << "vw - vk" << v_w - v_k << "\n"
                      << "vv - vk dot vw - vk: " << (v_v - v_k).dot(v_w - v_k) << "\n";
#endif
            double mod1 = mod((v_v - v_k).cross(v_w - v_k)),
                   mod2 = mod((v_v - v_l).cross(v_w - v_l));
            assert(mod1 != 0.0 && mod2 != 0.0);
            double cot_k = (v_v - v_k).dot(v_w - v_k) / mod1;
            double cot_l = (v_v - v_l).dot(v_w - v_l) / mod2;

            out << "mod1: " << mod1 << "\n";
            out << "mod2: " << mod2 << "\n";

            if ((mod((v_v - v_k).cross(v_w - v_k))) == 0.0) {
                // std::cout << "v: " << v_v << "\n";
                // std::cout << "k: " << v_k << "\n";
                std::cout << "------------v: " << v_v << "\n";
                std::cout << "------------w: " << v_w << "\n";
                std::cout << "------------k: " << v_k << "\n";
                std::cout << "------------l: " << v_l << "\n";
                // std::cout << "w - k: " << v_w - v_k << "\n";   
            }
            if ((mod((v_v - v_l).cross(v_w - v_l))) == 0.0) {
                // std::cout << "v: " << v_v << "\n";
                // std::cout << "k: " << v_l << "\n";
                std::cout << "------------v: " << v_v << "\n";
                std::cout << "------------w: " << v_w << "\n";
                std::cout << "------------k: " << v_k << "\n";
                std::cout << "------------l: " << v_l << "\n";     
            }
#if MESH_DEBUG
            std::cout << "cot_k: " << cot_k << "\n";
            std::cout << "cot_l: " << cot_l << "\n";
#endif
            return cot_k + cot_l;
        }
    }

    int num_of_boundary(CMyMesh & mesh)
    {
        int sum = 0;
        for (CMyMesh::MeshVertexIterator viter(&mesh); !viter.end(); ++viter)
        {
            CMyVertex *v = *viter;
            if (v->boundary()) {
                ++sum;
            }
        }
        return sum;
    }

    int num_of_inner(CMyMesh & mesh)
    {
        return mesh.vertices().size() - num_of_boundary(mesh);
    }

    double length_halfedge(CHalfEdge * he)
    {
        CPoint p1 = he->vertex()->point(),
               p2 = he->source()->point();
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
               (p1[1] - p2[1]) * (p1[1] - p2[1]) +
               (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }

    /**
 * [calculate the sum of hamonic energy of mesh]
 * @return [the sum of hamonic energy of mesh]
 */
    double harmonicEnergy(CMyMesh & mesh)               
    {

        double E = 0.0;
        for (std::map<CMyVertex*, Eigen::Vector3d>::iterator it = f.begin(); it != f.end(); ++ it) {
            CMyVertex * v = it->first;
            for (CMyMesh::VertexVertexIterator vviter(v); !vviter.end(); ++vviter) {
                CMyVertex * w = *vviter;
                double E_f = ((f[v][0] - f[w][0]) * (f[v][0] - f[w][0]) + (f[v][1] - f[w][1]) * (f[v][1] - f[w][1]));
                double k_v_w = k(v, w, mesh);
                // std::cout << "E_f: " << E_f << std::endl;
                //std::cout << "k: " << k_v_w << std::endl;
                E += k_v_w * E_f;
            }
        }
        return E;
    }

    void set_boundary_disk(CMyMesh & mesh)
    {
        /*some check, wheher genus is 1 or not*/
        int nv = mesh.numVertices();
        int ne = mesh.numEdges();
        int nf = mesh.numFaces();
        int euler_char = nv - ne + nf;
        CMyMesh::CBoundary boundary(&mesh);
        std::vector<CMyMesh::CLoop*> & loops = boundary.loops();
        int nb = loops.size();
        int genus = (2 - (euler_char + nb)) / 2;
        assert(genus == 0);
        /*end check*/

        /*
                init f
         */
        for (CMyMesh::MeshVertexIterator viter(&mesh); !viter.end(); ++viter) {
            CMyVertex * v = *viter;
            CPoint p = v->point();
            f.insert(std::make_pair(v, Eigen::Vector3d(0, 0, 0)));
            org.insert(std::make_pair(v, Eigen::Vector3d(p[0], p[1], p[2])));
        }

        std::list<CMyHalfEdge*> & boundary_list = loops[0]->halfedges();
        std::cout << "loops: " << boundary_list.size() << "\n";
        std::vector<CMyVertex* > _vlist_;   //boundary vertexs

        std::vector<double> vlength;   //[0, l01, l12, l23, l34,....]
        double total_length = 0;

        int i = 0;
        for (std::list<CMyHalfEdge*>::iterator it = boundary_list.begin(); it != boundary_list.end(); ++ it) {
            CMyHalfEdge *he = *it;
            _vlist_.push_back((CMyVertex*)he->vertex());
            if (vlength.size() == 0) {
                vlength.push_back(0.0);
            }
            else {
                vlength.push_back(length_halfedge(he) + vlength[i - 1]);
            }
            ++i;
            total_length += length_halfedge(he);
        }

        for(int i = _vlist_.size() - 1; i >= 0 ; --i) {
            double theta = vlength[i] / total_length * 2.0 * M_PI;
            Eigen::Vector3d g_v;
            // g_v << cos(theta), sin(theta), 0.0;
            g_v << sin(theta), cos(theta), 0.0; 
            g.insert(std::pair<CMyVertex*, Eigen::Vector3d>(_vlist_[i], g_v));
            f[_vlist_[i]] = g_v;
            _vlist_[i]->point() = CPoint(g_v[0], g_v[1], g_v[2]);
#if MESH_DEBUG
            std::cout << g_v << "\n";
#endif
        }

        std::cout << "energy: " << harmonicEnergy(mesh) << "\n";
    }

    void set_boundary_quad(CMyMesh & mesh)
    {
         /*some check, wheher genus is 1 or not*/
        int nv = mesh.numVertices();
        int ne = mesh.numEdges();
        int nf = mesh.numFaces();
        int euler_char = nv - ne + nf;
        CMyMesh::CBoundary boundary(&mesh);
        std::vector<CMyMesh::CLoop*> & loops = boundary.loops();
        int nb = loops.size();
        int genus = (2 - (euler_char + nb)) / 2;
        assert(genus == 0);
        /*end check*/

        std::list<CMyHalfEdge*> & boundary_list = loops[0]->halfedges();
        std::vector<CMyVertex* > _vlist_;
        std::vector<CMyVertex* > _vlist_1;          //(0, 0)~(1, 0)
        std::vector<CMyVertex* > _vlist_2;          //(1, 0)~(1, 1)
        std::vector<CMyVertex* > _vlist_3;          //(1, 1)~(0, 1) 
        std::vector<CMyVertex* > _vlist_4;          //(0, 1)~(0, 0)

        for (auto & he : boundary_list) {
            _vlist_.push_back((CMyVertex*)he->vertex());
        }

        /*
                init f
         */
        for (CMyMesh::MeshVertexIterator viter(&mesh); !viter.end(); ++viter) {
            CMyVertex * v = *viter;
            CPoint p = v->point();
            f.insert(std::make_pair(v, Eigen::Vector3d(0, 0, 0)));
            org.insert(std::make_pair(v, Eigen::Vector3d(p[0], p[1], p[2])));
        }

        int segmt = _vlist_.size() / 4;
//--------------------------------------------------

        for (int i = 0; i < segmt; ++i) {
            _vlist_1.push_back(_vlist_[i]);
        }

        for (int i = segmt; i < 2 * segmt; ++ i) {
            _vlist_2.push_back(_vlist_[i]);
        }

        for (int i = 2 * segmt; i < 3 * segmt; ++i) {
            _vlist_3.push_back(_vlist_[i]);
        }

        for (int i = 3 * segmt; i < _vlist_.size(); ++i) {
            _vlist_4.push_back(_vlist_[i]);
        }

        for(int i = 0; i < _vlist_1.size(); ++i) {          //  (0, 0) ~ (1, 0)
            Eigen::Vector3d v((double)i / (double)_vlist_1.size(), 0.0, 0.0);
            g.insert(std::make_pair(_vlist_1[i], v));
            f[_vlist_1[i]] = v;
            CPoint g_v(v[0], v[1], v[2]);
            _vlist_1[i]->point() = g_v;
#if MESH_DEBUG
            std::cout << g_v << "\n";
#endif
        }

        for(int i = 0; i < _vlist_2.size(); ++i) {          //  (0, 0) ~ (1, 0)
            Eigen::Vector3d v(1.0, (double)i / (double)_vlist_2.size(), 0.0);
            g.insert(std::make_pair(_vlist_2[i], v));
            f[_vlist_2[i]] = v;
            CPoint g_v(v[0], v[1], v[2]);
            _vlist_2[i]->point() = g_v;
#if MESH_DEBUG
            std::cout << g_v << "\n";
#endif
        }
        
        for(int i = 0; i < _vlist_3.size(); ++i) {          //  (0, 0) ~ (1, 0)
            Eigen::Vector3d v((double)(_vlist_3.size() - i) / (double)_vlist_3.size(), 1.0, 0.0);
            g.insert(std::make_pair(_vlist_3[i], v));
            f[_vlist_3[i]] = v;
            CPoint g_v(v[0], v[1], v[2]);
            _vlist_3[i]->point() = g_v;
#if MESH_DEBUG
            std::cout << g_v << "\n";
#endif
        }
        for(int i = 0; i < _vlist_4.size(); ++i) {          //  (0, 0) ~ (1, 0)
            Eigen::Vector3d v(0.0, (double)(_vlist_4.size() - i) / (double)_vlist_4.size(), 0.0);
            g.insert(std::make_pair(_vlist_4[i], v));
            f[_vlist_4[i]] = v;
            CPoint g_v(v[0], v[1], v[2]);
            _vlist_4[i]->point() = g_v;
#if MESH_DEBUG
            std::cout << g_v << "\n";
#endif
        }
    }

    void map(CMyMesh & mesh)                            /*--ALGORITHM:2--*/
    {
        // set_boundary(mesh);
        int num_boundary = num_of_boundary(mesh);
        int num_inner = num_of_inner(mesh);
#if MESH_DEBUG
        std::cout << "num_boundary: " << num_boundary << "\n"
                  << "num_inner: " << num_inner << "\n";
#endif

        // assert(num_boundary == g.size());        //all the functions of boundary vertex
        std::cout << "boundary: " << num_boundary << "\t" << "g_size: " << g.size() << "\t";

        double E = harmonicEnergy(mesh);        //init E    /*--##FORMULA:3--*/

        do {

            Energy = E;                                     /*--##FORMULA:5--*/
            for (std::map<CMyVertex*, Eigen::Vector3d>::iterator it = f.begin(); it != f.end(); ++ it) {
                CMyVertex * v = it->first;
                Eigen::Vector3d p = it->second;
                double sum_p_0 = 0, sum_p_1 = 0;
                double sum_k = 0;
                if (!v->boundary()) {
                    for (CMyMesh::VertexVertexIterator vviter(v); !vviter.end(); ++vviter) {
                        CMyVertex * w = *vviter;
                        double k_v_w = k(v, w, mesh); 
                        sum_p_0 += k_v_w * f[w][0];  // x
                        sum_p_1 += k_v_w * f[w][1];  // y
                        sum_k += k_v_w;
                    }
                    Eigen::Vector3d v_t(sum_p_0 / sum_k, sum_p_1 / sum_k, 0.0); /*--##FORMULA:7--*/
                    f[v] = v_t;
                    v->point() = CPoint(v_t[0], v_t[1], v_t[2]);
                }
            }
            E = harmonicEnergy(mesh);

            std::cout << "diff: " << fabs(Energy - E) << "\n"; 

        } while (fabs(Energy - E) > this->e);
    } 

};

void generateHarmornicMap(int type, CMyMesh & mesh)
{
    harmornicMap m(0.1);
    if (type == DISK) {
        m.set_boundary_disk(mesh);
        m.map(mesh);
    }
    else if (type == QUAD) {
        m.set_boundary_quad(mesh);
        m.map(mesh);
    }
}