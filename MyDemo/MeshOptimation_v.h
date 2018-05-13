#pragma once
#include "HarmonicMap_v.h"


namespace meshOptimation {

    using namespace MeshLib;

    double length_halfedge(CHalfEdge * he)
    {
        CPoint p1 = he->vertex()->point(),
               p2 = he->he_prev()->vertex()->point();
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
               (p1[1] - p2[1]) * (p1[1] - p2[1]) +
               (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }

    double miniFaceAngle(CFace * face) 
    {
        CHalfEdge *h0 = face->halfedge(),
                  *h1 = h0->he_next(),
                  *h2 = h0->he_prev();
        CPoint e_0_a = ((CMyVertex*)h1->vertex())->huv() - ((CMyVertex*)h0->vertex())->huv(),
               e_0_b = ((CMyVertex*)h2->vertex())->huv() - ((CMyVertex*)h0->vertex())->huv(),
               e_1_a = ((CMyVertex*)h0->vertex())->huv() - ((CMyVertex*)h1->vertex())->huv(),
               e_1_b = ((CMyVertex*)h2->vertex())->huv() - ((CMyVertex*)h1->vertex())->huv(),
               e_2_a = ((CMyVertex*)h0->vertex())->huv() - ((CMyVertex*)h2->vertex())->huv(),
               e_2_b = ((CMyVertex*)h1->vertex())->huv() - ((CMyVertex*)h2->vertex())->huv();
        double ang_0 = acos(e_0_a * e_0_b / mod(e_0_a) / mod(e_0_b)),
               ang_1 = acos(e_1_a * e_1_b / mod(e_1_a) / mod(e_1_b)),
               ang_2 = acos(e_2_a * e_2_b / mod(e_2_a) / mod(e_2_b));
        // std::cout << "angle0:" << ang_0 << "\t" << "angle1:" << ang_1 << "\t" << "angle2:" << ang_2 << "\n"; 
        return std::min(ang_0, std::min(ang_1, ang_2));
    }

    CPoint baryCenter(CFace * face)
    {
        CPoint  v0 = ((CMyVertex*)face->halfedge()->vertex())->huv(),
                v1 = ((CMyVertex*)face->halfedge()->he_next()->vertex())->huv(),
                v2 = ((CMyVertex*)face->halfedge()->he_prev()->vertex())->huv();
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

    void legalizeFaceByInsert(CMyMesh& mesh, CFace * theface) 
    {

        CPoint p = baryCenter(theface);    // calculate circumcircle-center

        CFace* face = LocatePoint(mesh, p);
        if (face != NULL) {                 // if face exist, insert the circumcircle-center
            // std::cout << "insert." << "\n";
            // printFaceInfo(face); 
            CHalfEdge *h0 = face->halfedge(),
                      *h1 = h0->he_next(),
                      *h2 = h0->he_prev();
            CPoint p_org_avg =  (h0->vertex()->point() + 
                                 h1->vertex()->point() + 
                                 h2->vertex()->point()) / 3.0;
            CPoint n_avg = (h0->vertex()->normal() + 
                                h1->vertex()->normal() + 
                                h2->vertex()->normal()) / 3.0;
            CVertex* _v_ = faceSplit(mesh, face, p, p_org_avg, n_avg);

            assert(h0->he_next()->vertex() == h1->he_next()->vertex() 
                   && h0->he_next()->vertex() == h2->he_next()->vertex());
            CVertex* v = h0->he_next()->vertex();
            legalizeEdge(mesh, v, h0->edge());
            legalizeEdge(mesh, v, h1->edge());
            legalizeEdge(mesh, v, h2->edge());    
        }
        else {                               // else, split the longest edge
            // std::cout << "out\n";
        }

    } 

    void legalizeFaceBySplit(CMyMesh & mesh, CFace * theface, double length_ref)
    {
        // std::cout << "begin spliting...\n";
        CHalfEdge *h0 = theface->halfedge(),
                *h1 = h0->he_next(),
                *h2 = h0->he_prev();
        double l0 = length_halfedge(h0),
                l1 = length_halfedge(h1),
                l2 = length_halfedge(h2);
        if (std::min(l2, std::min(l0, l1)) > length_ref / 2.0) {
            if (l0 > l1 && l0 > l2) {
                myEdgeSplit(mesh, h0->edge());
            }
            if (l1 > l0 && l1 > l2) {
                myEdgeSplit(mesh, h1->edge());
            }
            if (l2 > l1 && l2 > l0) {
                myEdgeSplit(mesh, h2->edge());
            }
        }
        // std::cout << "done split...\n";
    }

    void legalizeMesh(CMyMesh& mesh)  
    {
        std::list<CMyFace*>::iterator it = (mesh.faces()).begin();
        CFace * pF = *it;
        double min_length = ( length_halfedge(pF->halfedge()) + 
                            length_halfedge(pF->halfedge()->he_prev()) +
                            length_halfedge(pF->halfedge()->he_next()) ) / 3.0;
        int i = 0, j = 0;
        for (; i < mesh.faces().size(); ++i)
        {
            std::cout << "iterations: " << "\t" << i << "\n";
            std::list<CMyFace*>::iterator it = (mesh.faces()).begin(); 
            std::advance(it, i); 
            CFace * pF = *it;
            if (miniFaceAngle(pF) < 0.6) {    //30
                legalizeFaceByInsert(mesh, pF);
            }
            else {
                // std::cout << "pass." << "\n";  
            }
        }
        for (; j < mesh.faces().size(); ++j)
        {
            std::cout << "iterations: " << "\t" << i + j << "\n";
            std::list<CMyFace*>::iterator it = (mesh.faces()).begin(); 
            std::advance(it, j); 
            CFace * pF = *it;
            if (miniFaceAngle(pF) < 0.6) {    //30
                legalizeFaceBySplit(mesh, pF, min_length);
            }
            else {
                // std::cout << "pass." << "\n";
            }
        }

        // for (CMyMesh::MeshVertexIterator viter(&mesh); !viter.end(); ++viter)
        // {
        //     CMyVertex * pV = *viter;
        //     pV->point() = pV->huv();
        // } 
    }

}

