/*!
*      \file HarmonicMapper.h
*      \brief Algorithm for harmonic mapping
*	   \author David Gu
*      \date Document 10/07/2010
*
*		Simple harmonic map, that maps a topological disk to the unit disk 
*		with minimal harmonic energy.
*/

// Simple harmonic map, that maps topological disk to a unit disk using harmonic map.

#ifndef _HARMONIC_MAPPER_H_
#define _HARMONIC_MAPPER_H_

#include <vector>
#include "HarmonicMapperMesh.h"

#ifndef PI
#define PI 3.141592653589793238462643383279
#endif

namespace MeshLib
{
/*!
 *	\brief CHarmonicMapper class
 *
 *  Compute the harmonic map by solving Dirichlet problem
 * \f[
 *		\left\{
 *		\begin{array}{ccc}
 *		 \Delta u &\equiv& 0\\
 *		 u|_{\partial \Omega} &=& f
 *		 \end{array}
 *		\right.
 * \f]
 */
	class CHarmonicMapper
	{
	public:
		/*!	CHarmonicMapper constructor
		 *	\param pMesh the input mesh
		 */
		CHarmonicMapper(CHMMesh* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CHarmonicMapper();
		/*!  Compute the harmonic map using direct method
		 */
		void _map();
		/*!	Iterative method compute harmonic map
		 *	\param epsilon error threshould
		 */	
		void _iterative_map( double threshould = 5e-4 );

	protected:
		/*!	fix the boundary vertices to the unit circle
		 *  using arc length parameter
		 */
		void _set_boundary();
		/*!	The input surface mesh
		 */
		CHMMesh* m_pMesh;
		/*!	The boundary of m_pMesh
		 */
		CHMMesh::CBoundary m_boundary;
		
		/*!	number of interior vertices
		 */ 
		int m_interior_vertices;
		/*! number of boundary vertices
		*/
		int m_boundary_vertices;
	};
}

#endif

/*-------------------------------------------------------------------------------------------
//Example using CHarmonicMapper

#include "HarmonicMapper/HarmonicMapper.h"

using namespace MeshLib;	


void help(char * exe )
{
	printf("Usage:\n");
	printf("%s -harmonic_map input_mesh output_mesh\n", exe);


}
int main( int argc, char * argv[] )
{
	if( argc < 3 )
	{
		help( argv[0] );
		return 0;
	}
	if( strcmp( argv[1] , "-harmonic_map") == 0 )
	{
		CHMMesh mesh;
		mesh.read_m( argv[2] );

		CHarmonicMapper mapper( & mesh );
		mapper._map();
		mesh.write_m( argv[3] );
		return 0;
	}

	help( argv[0] );
	return 0;
} 
-------------------------------------------------------------------------------------------*/