#pragma once
#include "HarmonicMap.h"


namespace meshOptimation {

    using namespace MeshLib;
    harmornicMap m(0.3);

    void insertVertex(CMyMesh& mesh, CPoint p)
    {
        CFace* face = delaunayTriangulation::LocatePoint(mesh, p);
        if (face != NULL) {
            // printFaceInfo(face);
            CHalfEdge *h0 = face->halfedge(),
                      *h1 = h0->he_next(),
                      *h2 = h0->he_prev();
            delaunayTriangulation::faceSplit(mesh, face, p);
            assert(h0->he_next()->vertex() == h1->he_next()->vertex() 
                   && h0->he_next()->vertex() ==h2->he_next()->vertex());
            CVertex* v = h0->he_next()->vertex();
            delaunayTriangulation::legalizeEdge(mesh, v, h0->edge());
            delaunayTriangulation::legalizeEdge(mesh, v, h1->edge());
            delaunayTriangulation::legalizeEdge(mesh, v, h2->edge());                //problem!!!!!!!!!!!!
        } else {
            std::cout << p << "---out---\n";
        }
    }

    double miniFaceAngle(CFace * face) 
    {
        CHalfEdge *h0 = face->halfedge(),
                  *h1 = h0->he_next(),
                  *h2 = h0->he_prev();
        CPoint e_0_a = h1->vertex()->point() - h0->vertex()->point(),
               e_0_b = h2->vertex()->point() - h0->vertex()->point(),
               e_1_a = h0->vertex()->point() - h1->vertex()->point(),
               e_1_b = h2->vertex()->point() - h1->vertex()->point(),
               e_2_a = h0->vertex()->point() - h2->vertex()->point(),
               e_2_b = h1->vertex()->point() - h2->vertex()->point();
        double ang_0 = acos(e_0_a * e_0_b / mod(e_0_a) / mod(e_0_b)),
               ang_1 = acos(e_1_a * e_1_b / mod(e_1_a) / mod(e_1_b)),
               ang_2 = acos(e_2_a * e_2_b / mod(e_2_a) / mod(e_2_b));
        return std::min(ang_0, std::min(ang_1, ang_2));
    }

    CPoint baryCenter(CFace * face)
    {
        CPoint  v0 = face->halfedge()->vertex()->point(),
                v1 = face->halfedge()->he_next()->vertex()->point(),
                v2 = face->halfedge()->he_prev()->vertex()->point();
        // std::cout << "face: " << v0 << ", " << v1 << ", " << v2 << "\n";
        CPoint  va = (v0 + v1) / 2.0,
                vb = (v0 + v2) / 2.0;

        double ka = (v0[0] - v1[0]) / (v1[1] - v0[1]),
               kb = (v0[0] - v2[0]) / (v2[1] - v0[1]);
        double ba = va[1] - ka * va[0],
               bb = vb[1] - kb * vb[0];

        double x = (bb - ba) / (ka - kb),
               y = ka * (bb - ba) / (ka - kb) + ba;
        return CPoint(x, y, 0.0);
    }

    void legalizeFace(CMyMesh& mesh, CFace * theface) 
    {

        std::map<CMyVertex*, Eigen::Vector3d> & f = m.f;    //all the result 

        std::map<CMyVertex*, Eigen::Vector3d> & org = m.org;  //origin

        CPoint p = baryCenter(theface);
        // std::cout << "barycenter: " << p << "\n";

        CFace* face = delaunayTriangulation::LocatePoint(mesh, p);
        if (face != NULL) {
            std::cout << "insert." << "\n";
            // printFaceInfo(face); 
            CHalfEdge *h0 = face->halfedge(),
                      *h1 = h0->he_next(),
                      *h2 = h0->he_prev();
            Eigen::Vector3d v_org = (org[(CMyVertex*)h0->vertex()] + 
                                     org[(CMyVertex*)h1->vertex()] + 
                                     org[(CMyVertex*)h2->vertex()]) / 3.0;
            CVertex* _v_ = delaunayTriangulation::faceSplit(mesh, face, p);
            f.insert(std::make_pair<CMyVertex*, Eigen::Vector3d>((CMyVertex*)_v_, 
                                                                  Eigen::Vector3d(p[0], p[1], p[2])));
            org.insert(std::pair<CMyVertex*, Eigen::Vector3d>((CMyVertex*)_v_, v_org));
            assert(h0->he_next()->vertex() == h1->he_next()->vertex() 
                   && h0->he_next()->vertex() ==h2->he_next()->vertex());
            CVertex* v = h0->he_next()->vertex();
            delaunayTriangulation::legalizeEdge(mesh, v, h0->edge());
            delaunayTriangulation::legalizeEdge(mesh, v, h1->edge());
            delaunayTriangulation::legalizeEdge(mesh, v, h2->edge());                //problem!!!!!!!!!!!!
        }
        else {
            std::cout << "outof" << "\n";
        }

    }

    void legailizeMesh(CMyMesh& mesh)
    {
        
        // m.set_boundary_disk(mesh);
        // m.map(mesh);

        for (int i = 0; i < mesh.faces().size(); ++i)
        {
            std::cout << "iterations: " << "\t" << i << "\n";
            std::list<CMyFace*>::iterator it = (mesh.faces()).begin(); 
            std::advance(it, i);
            CFace * pF = *it;
            if (miniFaceAngle(pF) < 0.6) {    //30
                legalizeFace(mesh, pF);
            }
            else {
                std::cout << "pass." << "\n";
            }
        }

        // for (CMyMesh::MeshVertexIterator viter(&mesh); !viter.end(); ++viter)
        // {
        //     CVertex * pV = *viter;
        //     Eigen::Vector3d v = m.org[(CMyVertex*)pV];
        //     pV->point() = CPoint(v[0], v[1], v[2]);
        // }
        // for (int i = 0; i < mesh.faces().size(); ++i)
        // {
        //     CFace * pF = mesh.faces()[i];
        //     if (miniFaceAngle(pF) < 0.5233) {    //30
        //         legalizeFace(pF);
        //     }
        // }
    }

    CVertex* faceSplit(CMyMesh &mesh, CFace* f, CPoint p, CPoint p_org, CPoint n_avg)
    {
        int ivid = mesh.vertices().size() + 1;
        CVertex *v = mesh.createVertex(ivid);
        ((CMyVertex*)v)->huv() = p;
        v->point() = p_org;
        v->normal() = n_avg;
        v->boundary() = false;

        int ifid = mesh.faces().size() + 1;
        std::vector<CHalfEdge*> hes0, hes1, hes2;
        CHalfEdge *he0 = f->halfedge(),
                  *he1 = he0->he_next(),
                  *he2 = he0->he_prev();

//-----------------------------------------------f0 = f
        CHalfEdge *hes0_1 = new CHalfEdge(),
                  *hes0_2 = new CHalfEdge();
        hes0_1->vertex() = v;
        v->halfedge() = hes0_1;
        hes0_2->vertex() = he2->vertex();
        he2->vertex()->halfedge() = hes0_2;

        hes0.push_back(he0);
        hes0.push_back(hes0_1);
        hes0.push_back(hes0_2);
        //linking to each other
        for(int i = 0; i < 3; i ++ )
        {
            hes0[i]->he_next() = hes0[(i+1)%3];
            hes0[i]->he_prev() = hes0[(i+2)%3];
        }
        //connection with edge
        CEdge *e1 = mesh.createEdge((CMyVertex*)he0->vertex(), (CMyVertex*)v),
              *e2 = mesh.createEdge((CMyVertex*)v, (CMyVertex*)he2->vertex());
        e1->halfedge(0) = hes0_1;
        e2->halfedge(0) = hes0_2;
        hes0_1->edge() = e1;
        hes0_2->edge() = e2;
        //linking to face
        for(int i = 0; i < 3; i ++ )
        {
            hes0[i]->face()   = f;
            f->halfedge()    = hes0[i];
        }
//-----------------------------------------------f1
        CMyFace* f1 = new CMyFace();
        assert (f1 != NULL);
        f1->id() = ifid;
        mesh.faces().push_back(f1);
        mesh.map_face().insert( std::pair<int, CMyFace*>(ifid++, f1) );

        CHalfEdge *hes1_1 = new CHalfEdge(),
                  *hes1_2 = new CHalfEdge();
        hes1_1->vertex() = v;
        v->halfedge() = hes1_1;
        hes1_2->vertex() = he0->vertex();
        he0->vertex()->halfedge() = hes1_2;        //problem???

        hes1.push_back(he1);
        hes1.push_back(hes1_1);
        hes1.push_back(hes1_2);
        //linking to each other
        for(int i = 0; i < 3; i ++ )
        {
            hes1[i]->he_next() = hes1[(i+1)%3];
            hes1[i]->he_prev() = hes1[(i+2)%3];
        }
        //connection with edge
        e1 = mesh.createEdge((CMyVertex*)he1->vertex(), (CMyVertex*)v);
        e2 = mesh.createEdge((CMyVertex*)v, (CMyVertex*)he0->vertex());
        assert (e1->halfedge(0) == NULL);
        e1->halfedge(0) = hes1_1;
        e2->halfedge(1) = hes1_2;
        hes1_1->edge() = e1;
        hes1_2->edge() = e2;
        //linking to face
        for(int i = 0; i < 3; i ++ )
        {
            hes1[i]->face()   = f1;
            f1->halfedge()    = hes1[i];
        }
//-----------------------------------------------f2
        CMyFace* f2 = new CMyFace();
        assert (f2 != NULL);
        f2->id() = ifid;
        mesh.faces().push_back(f2);
        mesh.map_face().insert( std::pair<int, CMyFace*>(ifid, f2) );

        CHalfEdge *hes2_1 = new CHalfEdge(),
                  *hes2_2 = new CHalfEdge();
        hes2_1->vertex() = v;
        v->halfedge() = hes2_1;
        hes2_2->vertex() = he1->vertex();
        he1->vertex()->halfedge() = hes2_2;        //problem???

        hes2.push_back(he2);
        hes2.push_back(hes2_1);
        hes2.push_back(hes2_2);
        //linking to each other
        for(int i = 0; i < 3; i ++ )
        {
            hes2[i]->he_next() = hes2[(i+1)%3];
            hes2[i]->he_prev() = hes2[(i+2)%3];
        }
        //connection with edge
        e1 = mesh.createEdge((CMyVertex*)he2->vertex(), (CMyVertex*)v);
        e2 = mesh.createEdge((CMyVertex*)v, (CMyVertex*)he1->vertex());
        assert (e1->halfedge(1) == NULL);
        e1->halfedge(1) = hes2_1;
        e2->halfedge(1) = hes2_2;
        hes2_1->edge() = e1;
        hes2_2->edge() = e2;
        //linking to face
        for(int i = 0; i < 3; i ++ )
        {
            hes2[i]->face()   = f2;
            f2->halfedge()    = hes2[i];
        }

        return v;
    }

   void edgeSwap(CMyEdge* e)          //not boundary
    {
        assert(e->boundary() == false);
        CHalfEdge *oldhe1 = e->halfedge(0),
                  *oldhe2 = e->halfedge(1);
        //-----------------------------------------------------------------10-07
        CVertex *oldv1 = oldhe1->vertex(),
                *oldv2 = oldhe2->vertex();
        CVertex* pV = ( oldv1->id() < oldv2->id() ) ? oldv1 : oldv2;
        std::list<CEdge*> & oldledges = pV->edges();
        std::list<CEdge*>::iterator it =  std::find(oldledges.begin(), oldledges.end(), e);
        if (it != oldledges.end()) {
            // std::cout << "erase!" << "\n";
            oldledges.erase(it);
        }
        //-----------------------------------------------------------------10-07
        

        CFace *f1 = oldhe1->face(),
              *f2 = oldhe2->face();

        std::vector<CHalfEdge*> hes1, hes2;
        hes1.push_back(oldhe1->he_prev());
        hes1.push_back(oldhe1->he_sym()->he_next());

        hes2.push_back(oldhe2->he_prev());
        hes2.push_back(oldhe2->he_sym()->he_next());

        CVertex *v1 = oldhe1->he_next()->vertex(),
                *v2 = oldhe2->he_next()->vertex();
        CHalfEdge *newhe1 = new CHalfEdge(),
                  *newhe2 = new CHalfEdge();
        newhe1->vertex() = v1;
        v1->halfedge() = newhe1;
        hes1.push_back(newhe1);
        newhe2->vertex() = v2;
        v2->halfedge() = newhe2;
        hes2.push_back(newhe2);

        //linking to each other
        for(int i = 0; i < 3; i ++ )
        {
            hes1[i]->he_next() = hes1[(i+1)%3];
            hes1[i]->he_prev() = hes1[(i+2)%3];
        }
        //linking to each other
        for(int i = 0; i < 3; i ++ )
        {
            hes2[i]->he_next() = hes2[(i+1)%3];
            hes2[i]->he_prev() = hes2[(i+2)%3];
        }

        //connection to edge
        CMyEdge *&olded = e;
        olded->halfedge(0) = newhe1;
        olded->halfedge(1) = newhe2;
        newhe1->edge() = olded;
        newhe2->edge() = olded;
        //------------------------------------------------------------------10-07
        CVertex* pV_ = ( v1->id() < v2->id() ) ? v1 : v2;
        std::list<CEdge*> & ledges = pV_->edges();
        ledges.push_back(olded);
        //------------------------------------------------------------------10-07

        //linking to face
        for(int i = 0; i < 3; i ++ )
        {
            hes1[i]->face()   = f1;
            f1->halfedge()    = hes1[i];
        }
        //linking to face
        for(int i = 0; i < 3; i ++ )
        {
            hes2[i]->face()   = f2;
            f2->halfedge()    = hes2[i];
        }
    }
}