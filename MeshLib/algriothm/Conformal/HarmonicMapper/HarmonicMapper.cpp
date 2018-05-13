/*!
*      \file HarmonicMapper.cpp
*      \brief Implement CHarmonicMapper class
*	   \author David Gu
*      \date Document 10/07/2010
*
*		Simple Harmonic Map, that maps topological disk to a unit disk using harmonic mapping.
*
*/

//#define _HARMONIC_MAP_DEBUG_

#include "HarmonicMapper.h"
#include "../Structure/Structure.h"
#include "../../../core/Mesh/iterators.h"
#include "Eigen/Sparse"

using namespace MeshLib;
using std::vector;

//Constructor
//Count the number of interior vertices
//Count the number of boundary vertices
//Compute edge weight

/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
CHarmonicMapper::CHarmonicMapper( CHMMesh* pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{

	int vid  = 0; //interior vertex ID 
	int vfid = 0; //boundary vertex ID

	for( CHMMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		CHarmonicVertex * pV = *viter;
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
	CStructure<CHMMesh, CHarmonicVertex, CHarmonicEdge, CHarmonicFace, CHarmonicHalfEdge> pC( m_pMesh );
	pC._embedding_2_metric();  //convert embedding to metric
	pC._metric_2_angle();	   //convert metric to angle
	pC._angle_2_Laplace();	   //convert angle to cotangent edge weight

#ifdef _HARMONIC_MAP_DEBUG_

	for( CHMMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		CHarmonicEdge * pE = *eiter;
		pE->weight() = 1.0;
	}
#endif
}

//Destructor
/*!
 *	CHarmonicMapper destructor
 */
CHarmonicMapper::~CHarmonicMapper()
{
}


//Set the boundary vertices to the unit circle
/*!
 *	Fix the boundary using arc length parameter
 */
void CHarmonicMapper::_set_boundary()
{
	//get the boundary half edge loop
	std::vector<CHMMesh::CLoop*> & pLs =  m_boundary.loops();
	assert( pLs.size() == 1 );
	CHMMesh::CLoop * pL = pLs[0];
	std::list<CHarmonicHalfEdge*> & pHs = pL->halfedges();
	
	//compute the total length of the boundary
	double sum = 0;
	for( std::list<CHarmonicHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
	{
		CHarmonicHalfEdge * ph = *hiter;
		CHarmonicEdge * pE = m_pMesh->halfedgeEdge( ph );
		sum += pE->length();
	}

	//parameterize the boundary using arc length parameter
	double l = 0;
	for( std::list<CHarmonicHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
	{
		CHarmonicHalfEdge * ph = *hiter;
		CHarmonicEdge * pE = m_pMesh->halfedgeEdge( ph );
		l += pE->length();
		double ang = l/sum * 2.0 * PI;
		CHarmonicVertex * pV = m_pMesh->halfedgeTarget( ph );
		pV->huv()= CPoint2( (cos(ang)+1.0)/2.0, (sin(ang)+1.0)/2.0 );
	}
}

//Compute the harmonic map with the boundary condition, direct method
/*!	Compute harmonic map using direct method
*/
void CHarmonicMapper::_map()
{
	//fix the boundary
	_set_boundary();

	std::vector<Eigen::Triplet<double> > A_coefficients;
	std::vector<Eigen::Triplet<double> > B_coefficients;

	
	//set the matrix A
	for( CHMMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		CHarmonicVertex * pV = *viter;
		if( pV->boundary() ) continue;
		int vid = pV->idx();

		double sw = 0;
		for( CHMMesh::VertexVertexIterator witer( pV ); !witer.end(); ++ witer )
		{
			CHarmonicVertex * pW = *witer;
			int wid = pW->idx();
			
			CHarmonicEdge * e = m_pMesh->vertexEdge( pV, pW );
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
		for( CHMMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		{
			CHarmonicVertex * pV = *viter;
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
		for( CHMMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		{
			CHarmonicVertex * pV = *viter;
			if( pV->boundary() ) continue;
			int id = pV->idx();
			pV->huv()[k] = x(id);
		}

	}
}

//Compute the harmonic map with the boundary condition, iterative method
/*!	Iterative method compute harmonic map
*	\param epsilon error threshould
*/
void CHarmonicMapper::_iterative_map( double epsilon )
{
	//fix the boundary
	_set_boundary();

	//move interior each vertex to its center of neighbors
	for( CHMMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		CHarmonicVertex * pV = *viter;
		if( pV->boundary() ) continue;
		
		pV->huv() = CPoint2(0,0);
	}
	
	while( true )
	{
		double error = -1e+10;
		//move interior each vertex to its center of neighbors
		for( CHMMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		{
			CHarmonicVertex * pV = *viter;
			if( pV->boundary() ) continue;
			
			double  sw = 0;
			CPoint2 suv(0,0);
			for( CHMMesh::VertexVertexIterator vviter(pV); !vviter.end(); vviter ++ )
			{
				CHarmonicVertex * pW = *vviter;
				CHarmonicEdge   * pE = m_pMesh->vertexEdge( pV, pW );
				double w = pE->weight();
				sw += w;
				suv = suv + pW->huv() * w;
			}
			suv /= sw;
			double verror = (pV->huv()-suv).norm();
			error = (verror > error )?verror:error; 
			pV->huv() = suv;
		}
#ifdef _HARMONIC_MAP_DEBUG_
		printf("Current max error is %f\n", error );
#endif
		if( error < epsilon ) break;
	}	
}

