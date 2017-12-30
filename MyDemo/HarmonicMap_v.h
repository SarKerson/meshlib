#pragma once
#include "MyMesh.h"
#include "../MeshLib/algriothm/Structure/Structure.h"
#include "../MeshLib/core/Mesh/iterators.h"
#include <eigen3/Eigen/Sparse>

using namespace MeshLib;

class harmornicMap
{
protected:
    /*! fix the boundary vertices to the unit circle
     *  using arc length parameter
     */
    void _set_boundary() 
    {
        //get the boundary half edge loop
        std::vector<CMyMesh::CLoop*> & pLs =  m_boundary.loops();
        // assert( pLs.size() == 1 );
        CMyMesh::CLoop * pL = pLs[0];
        std::list<CMyHalfEdge*> & pHs = pL->halfedges();
        
        //compute the total length of the boundary
        double sum = 0;
        for( std::list<CMyHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
        {
            CMyHalfEdge * ph = *hiter;
            CMyEdge * pE = m_pMesh->halfedgeEdge( ph );
            sum += pE->length();
        }

        //parameterize the boundary using arc length parameter
        double l = 0;
        for( std::list<CMyHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
        {
            CMyHalfEdge * ph = *hiter;
            CMyEdge * pE = m_pMesh->halfedgeEdge( ph );
            l += pE->length();
            double ang = l/sum * 2.0 * PI;
            CMyVertex * pV = m_pMesh->halfedgeTarget( ph );
            pV->huv()= CPoint( (cos(ang)+1.0)/2.0, (sin(ang)+1.0)/2.0, 0 );
        }
    }
    /*! The input surface mesh
     */
    CMyMesh* m_pMesh;
    /*! The boundary of m_pMesh
     */
    CMyMesh::CBoundary m_boundary;
    
    /*! number of interior vertices
     */ 
    int m_interior_vertices;
    /*! number of boundary vertices
    */
    int m_boundary_vertices;
public:
    /*! CHarmonicMapper constructor
     *  \param pMesh the input mesh
     */
    harmornicMap(CMyMesh* pMesh): m_pMesh(pMesh), m_boundary(m_pMesh)
    {
        int vid  = 0; //interior vertex ID 
        int vfid = 0; //boundary vertex ID

        for( CMyMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
        {
            CMyVertex * pV = *viter;
            if( pV->boundary() )
            {
                pV->idx() = vfid ++;
            }
            else
            {
                pV->idx() = vid  ++;
            }
        }

        m_interior_vertices = vid;
        m_boundary_vertices = vfid;
        
        //Compute cotangent edge weight
        CStructure<CMyMesh, CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> pC( m_pMesh );
        pC._embedding_2_metric();  //convert embedding to metric
        pC._metric_2_angle();      //convert metric to angle
        pC._angle_2_Laplace();     //convert angle to cotangent edge weight
    }
    /*! CHarmonicMapper destructor
     */
    ~harmornicMap() {}
    /*!  Compute the harmonic map using direct method
     */
    void _map() 
    {
        //fix the boundary
        _set_boundary();

        std::vector<Eigen::Triplet<double> > A_coefficients;
        std::vector<Eigen::Triplet<double> > B_coefficients;

        
        //set the matrix A
        for( CMyMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
        {
            CMyVertex * pV = *viter;
            if( pV->boundary() ) continue;
            int vid = pV->idx();

            double sw = 0;
            for( CMyMesh::VertexVertexIterator witer( pV ); !witer.end(); ++ witer )
            {
                CMyVertex * pW = *witer;
                int wid = pW->idx();
                
                CMyEdge * e = m_pMesh->vertexEdge( pV, pW );
                double w = e->weight();

                if( pW->boundary() )
                {
                    Eigen::Triplet<double> e(vid,wid,w);
                    B_coefficients.push_back( e );
                }
                else
                {
                    A_coefficients.push_back( Eigen::Triplet<double>(vid,wid, -w) );
                    //A_coefficients.push_back( Eigen::triplet<double> e(vid,wid,-w) );
                }
                sw += w;
            }
            //A.insert( vid, vid) = sw;
            A_coefficients.push_back( Eigen::Triplet<double>(vid,vid, sw ) );
        }


        Eigen::SparseMatrix<double> A( m_interior_vertices, m_interior_vertices );
        A.setZero();

        Eigen::SparseMatrix<double> B( m_interior_vertices, m_boundary_vertices );
        B.setZero();
        A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
        B.setFromTriplets(B_coefficients.begin(), B_coefficients.end());


        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
        std::cerr << "Eigen Decomposition" << std::endl;
        solver.compute(A);
        std::cerr << "Eigen Decomposition Finished" << std::endl;
        
        if( solver.info() != Eigen::Success )
        {
            std::cerr << "Waring: Eigen decomposition failed" << std::endl;
        }


        for( int k = 0; k < 2; k ++ )
        {
            Eigen::VectorXd b(m_boundary_vertices);
            //set boundary constraints b vector
            for( CMyMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
            {
                CMyVertex * pV = *viter;
                if( !pV->boundary() ) continue;
                int id = pV->idx();
                b(id) = pV->huv()[k];
            }

            Eigen::VectorXd c(m_interior_vertices);
            c = B * b;

            Eigen::VectorXd x = solver.solve(c);
            if( solver.info() != Eigen::Success )
            {
                std::cerr << "Waring: Eigen decomposition failed" << std::endl;
            }

            //set the images of the harmonic map to interior vertices
            for( CMyMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
            {
                CMyVertex * pV = *viter;
                if( pV->boundary() ) continue;
                int id = pV->idx();
                pV->huv()[k] = x(id);
                // pV->point()[k] = x(id);   //add ----------------------------
            }

        }

        for( CMyMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
        {
            CMyVertex * pV = *viter;
            if( pV->boundary() ) continue;
            pV->huv()[2] = 0;
            // pV->point()[k] = x(id);   //add ----------------------------
        }
        std::cout << "map done...\n";
    }
    /*! Iterative method compute harmonic map
     *  \param epsilon error threshould
     */ 
    void _iterative_map( double threshould = 5e-4 )
    {
        //fix the boundary
        _set_boundary();

        //move interior each vertex to its center of neighbors
        for( CMyMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
        {
            CMyVertex * pV = *viter;
            if( pV->boundary() ) continue;
            
            pV->huv() = CPoint(0,0,0);
        }
        
        int i = 0;
        while( true )
        {
            std::cout << "iteration: " << ++i << "\n";
            double error = -1e+10;
            //move interior each vertex to its center of neighbors
            for( CMyMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
            {
                CMyVertex * pV = *viter;
                if( pV->boundary() ) continue;
                
                double  sw = 0;
                CPoint suv(0,0,0);
                for( CMyMesh::VertexVertexIterator vviter(pV); !vviter.end(); vviter ++ )
                {
                    CMyVertex * pW = *vviter;
                    CMyEdge   * pE = m_pMesh->vertexEdge( pV, pW );
                    double w = pE->weight();
                    sw += w;
                    suv = suv + pW->huv() * w;
                }
                suv /= sw;
                double verror = (pV->huv()-suv).norm();
                error = (verror > error )?verror:error; 
                pV->huv() = suv;
            }
            if( error < threshould ) break;
        }
    }



};

void generateHarmornicMap(CMyMesh & mesh)
{
    harmornicMap mapper(&mesh);
    // mapper._map();
    mapper._iterative_map();
}