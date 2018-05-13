/*!
*      \file GeneralHarmonicMapper.cpp
*      \brief Implement CGeneralHarmonicMapper class
*	   \author David Gu
*      \date Document 11/19/2012
*
*		Simple Harmonic Map, that maps topological disk to a unit disk using harmonic mapping.
*
*/


#include "GeneralHarmonicMapper.h"
#include "../Structure/Structure.h"
#include "../../../core/Mesh/iterators.h"

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
CGeneralHarmonicMapper::CGeneralHarmonicMapper( CHMMesh* pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
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

}

//Destructor
/*!
 *	CGeneralHarmonicMapper destructor
 */
CGeneralHarmonicMapper::~CGeneralHarmonicMapper()
{
}


//Set the boundary vertices to the unit circle
/*!
 *	Fix the boundary using arc length parameter
 */
void CGeneralHarmonicMapper::_set_boundary()
{
	//get the boundary half edge loop
	std::vector<CHMMesh::CLoop*> & pLs =  m_boundary.loops();
	assert( pLs.size() == 1 );
	CHMMesh::CLoop * pL = pLs[0];
	std::list<CHarmonicHalfEdge*> & pHs = pL->halfedges();

	std::list<CHarmonicHalfEdge*>  hedges;

	std::list<CHarmonicVertex*> corners;
	for( std::list<CHarmonicHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
	{
		CHarmonicHalfEdge * ph = *hiter;
		CHarmonicVertex   * pv = m_pMesh->halfedgeTarget( ph );
		
		int valence = 0;
		for( CHMMesh::VertexEdgeIterator veiter( pv ); !veiter.end(); veiter ++ )
		{
			CHarmonicEdge * pE = *veiter;
			valence += ( pE->sharp() )?1:0;
		}

		if( valence > 2 ) corners.push_back( pv );

		hedges.push_back( ph );
	}

	printf("Number of corners is %d\n", corners.size());
	
	int k = 1;
	for( std::list<CHarmonicVertex*>::iterator viter = corners.begin(); viter != corners.end(); viter ++ )
	{
		CHarmonicVertex * pv = *viter;
		CPoint2 p = CPoint2( pv->point()[0], pv->point()[1] );
		double r = 1.0/(k+1.0);
		p = p * ( 1 - r );
		pv->huv() = p;
		k++;
	}


	while( true )
	{
		CHarmonicHalfEdge * ph = hedges.front();
		CHarmonicVertex   * pv = m_pMesh->halfedgeSource( ph );
		std::list<CHarmonicVertex*>::iterator pos = std::find( corners.begin(), corners.end(), pv );
		if( pos != corners.end() ) break;

		hedges.pop_front();
		hedges.push_back( ph );
	}	

	std::vector<std::list<CHarmonicHalfEdge*>*> boundaries;
	std::list<CHarmonicHalfEdge*>* pSeg = NULL;

	while( !hedges.empty() )
	{
		CHarmonicHalfEdge * ph = hedges.front();
		hedges.pop_front();

		CHarmonicVertex * pv = m_pMesh->halfedgeSource( ph );
		std::list<CHarmonicVertex*>::iterator pos = std::find( corners.begin(), corners.end(), pv );
		if( pos != corners.end() )
		{
			pSeg = new std::list<CHarmonicHalfEdge*>;
			boundaries.push_back( pSeg );		
		}
		pSeg->push_back( ph );
	}

	printf("Number of segments is %d\n", boundaries.size());

	
	for( size_t i = 0; i < boundaries.size(); i ++ )
	{
		std::list<CHarmonicHalfEdge*>* pL = boundaries[i];

		double total_length = 0;
		for( std::list<CHarmonicHalfEdge*>::iterator hiter = pL->begin(); hiter != pL->end(); hiter ++ )
		{
			CHarmonicHalfEdge* ph = *hiter;
			CHarmonicEdge    * pe = m_pMesh->halfedgeEdge( ph );
			total_length += pe->length();
		}
		
		CHarmonicVertex * pSource = m_pMesh->halfedgeSource( pL->front() );
		CHarmonicVertex * pTarget = m_pMesh->halfedgeTarget( pL->back()  );

		double partial_sum = 0;
		for( std::list<CHarmonicHalfEdge*>::iterator hiter = pL->begin(); hiter != pL->end(); hiter ++ )
		{
			CHarmonicHalfEdge* ph = *hiter;
			CHarmonicEdge    * pe = m_pMesh->halfedgeEdge( ph );
			CHarmonicVertex   *pv = m_pMesh->halfedgeTarget(ph);
			partial_sum += pe->length();

			double ratio = partial_sum/total_length;
			pv->huv() = pSource->huv() * (1-ratio) + pTarget->huv() * ratio;
		}
	}
}

//Compute the harmonic map with the boundary condition, direct method
/*!	Compute harmonic map using direct method
*/
void CGeneralHarmonicMapper::_map()
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
				B_coefficients.push_back( Eigen::Triplet<double>(vid,wid,w)); 
			}
			else
			{
				A_coefficients.push_back( Eigen::Triplet<double>(vid,wid,-w)); 
			}
			sw += w;
		}

		A_coefficients.push_back( Eigen::Triplet<double>(vid,vid,sw)); 
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
		Eigen::VectorXd b( m_boundary_vertices );

		//set boundary constraints b vector
		for( CHMMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		{
			CHarmonicVertex * pV = *viter;
			if( !pV->boundary() ) continue;
			int id = pV->idx();
			b(id)  = pV->huv()[k];
		}
		
		Eigen::VectorXd c = B * b;

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

