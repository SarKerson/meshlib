/*!
*      \file GeneralHarmonicMapper.h
*      \brief Algorithm for general harmonic mapping
*	   \author David Gu
*      \date Document 11/19/2012
*
*		Simple harmonic map, that maps a topological disk to the unit disk 
*		with minimal harmonic energy.
*/

// Simple harmonic map, that maps topological disk to a unit disk using harmonic map.

#ifndef _GENERAL_HARMONIC_MAPPER_H_
#define _GENERAL_HARMONIC_MAPPER_H_

#include <vector>
#include "HarmonicMapperMesh.h"
#include <Eigen/Sparse>

#ifndef PI
#define PI 3.141592653589793238462643383279
#endif

namespace MeshLib
{
/*!
 *	\brief CGeneralHarmonicMapper class
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
	class CGeneralHarmonicMapper
	{
	public:
		/*!	CGeneralHarmonicMapper constructor
		 *	\param pMesh the input mesh
		 */
		CGeneralHarmonicMapper(CHMMesh* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CGeneralHarmonicMapper();
		/*!  Compute the harmonic map using direct method
		 */
		void _map();

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

