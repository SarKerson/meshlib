# meshlib
a halfedge-based mesh library for 3d-mesh 

You can find the core library in the **/MeshLib**, and demo in the **/MyDemo**; some mesh data for testing is in the **/Data**.
**/MyDemo/main.cpp** shows how to use this library.
in **/MyDemo/MyMesh.h**, there are some namespaces that implement some algrithoms.

# SOME TRIPS FOR USING
## WHAT SHOULD BE DONE WHILE BUILDING A MESH
### 1. label boundary for edges, halfedges, vertex
	//Label boundary edges
	for( typename std::list<CEdge*>::iterator eiter= m_edges.begin() ; eiter != m_edges.end() ; ++ eiter )
	{
		CEdge *     edge = *eiter;
		CHalfEdge * he[2];

		he[0] = edgeHalfedge( edge, 0 );
		he[1] = edgeHalfedge( edge, 1 );
		
		assert( he[0] != NULL );
		

		if( he[1] != NULL )
		{
			assert( he[0]->target() == he[1]->source() && he[0]->source() == he[1]->target() );

			if( he[0]->target()->id() < he[0]->source()->id() )
			{
				edge->halfedge(0 ) = he[1];
				edge->halfedge(1 ) = he[0];
			}

			assert( edgeVertex1(edge)->id() < edgeVertex2(edge)->id() );
		}
		else
		{
			he[0]->vertex()->boundary() = true;
			he[0]->he_prev()->vertex()->boundary()  = true;
		}

	}
    
### 2.remove the dangling vertexs
    typename std::list<CVertex*> dangling_vert;
	for(typename std::list<CVertex*>::iterator viter = m_verts.begin();  viter != m_verts.end() ; ++ viter )
	{
		CVertex *     v = *viter;
		if( v->halfedge() != NULL ) continue;
		dangling_verts.push_back( v );
	}

	for( typename std::list<CVertex*>::iterator  viter = dangling_verts.begin() ; viter != dangling_verts.end(); ++ viter )
	{
		CVertex * v = *viter;
		m_verts.remove( v );
		delete v;
		v = NULL;
	}
### 3.Arrange the boundary half_edge
  //Arrange the boundary half_edge of boundary vertices, to make its halfedge
	//to be the most ccw in half_edge

	for(typename std::list<CVertex*>::iterator viter = m_verts.begin();  viter != m_verts.end() ; ++ viter )
	{
		CVertex *     v = *viter;
		if( !v->boundary() ) continue;

		CHalfEdge * he = vertexMostCcwInHalfEdge( v );
		while( he->he_sym() != NULL )
		{
			he =  vertexNextCcwInHalfEdge ( he );
		}
		v->halfedge() = he;
	}

## WHAT SHOULD BE DONE WHILE CREATING NEW FACES
### 1. create new face and push it into m_faces, insert a pair <id, face> into m_map_face
    CFace * f = new CFace();
    assert( f != NULL );
    f->id() = id;
    m_faces.push_back( f );
    m_map_face.insert( std::pair<int,tFace>(id,f) );
    
### 2. create new halfedges of faces, binding with vertex
    //create halfedges
    std::vector<tHalfEdge> hes;
    for(size_t i = 0; i < v.size(); i ++ )
    {
        tHalfEdge pH = new CHalfEdge;
        assert( pH );
        CVertex * vert =  v[i];
        pH->vertex() = vert;
        vert->halfedge() = pH;
        hes.push_back( pH );
    }
### 3. linking halfeges to each other
    //linking to each other
    for( size_t i = 0; i < hes.size(); i ++ )
    {
        hes[i]->he_next() = hes[(i+1)%hes.size()];
        hes[i]->he_prev() = hes[(i+hes.size()-1)%hes.size()];
    }
### 4. linking halfedges to the face
    //linking to face
    for(size_t i = 0; i < hes.size(); i ++ )
    {
        hes[i]->face()   = f;
        f->halfedge()    = hes[i];
    }
### 5. linking halfedge to the edge
    //connecting with edge
    for( size_t i = 0; i < hes.size(); i ++ )
    {
        tEdge e = createEdge( v[i], v[(i+hes.size()-1)%hes.size()] );
        if( e->halfedge(0)  == NULL )
        {
            e->halfedge(0) = hes[i];
        }
        else
        {
            assert( e->halfedge(1) == NULL );
            if( e->halfedge(1) != NULL )
            {
                std::cout << "Illegal Face Construction " << id << std::endl;
            }
            e->halfedge(1) = hes[i];
        }
        hes[i]->edge() = e;
    }

## WHAT SHOULD BE DONE WHILE CHANGING EDGE
Because the call of the method BaseMesh.h : createEdge( tVertex  v1, tVertex  v2 )
edges of the vertex would be searched first, if there exit a edge that link to [v1, v2], then return it.

  tVertex pV = ( v1->id()<v2->id())?v1:v2;
	typename std::list<CEdge*> & ledges = (typename std::list<CEdge*> &) pV->edges();

	
	for(typename std::list<CEdge*>::iterator te = ledges.begin(); te != ledges.end(); te ++ )
	{
		CEdge	  * pE = *te;
		CHalfEdge * pH = (CHalfEdge*) pE->halfedge(0);
	
		if( pH->source() == v1 && pH->target() == v2 ) 
		{
			return pE;		
		}
		if( pH->source() == v2 && pH->target() == v1 )
		{
			return pE;
		}
	}

else, create a new edge and push it into the list of pV(v1 or v2)

	//new edge
	CEdge * e = new CEdge;
	assert( e != NULL );
	m_edges.push_back( e );
	e->id() = (int)m_edges.size();
	ledges.push_back( e );
### removing a edge

    CVertex *oldv1 = oldhe1->vertex(),
            *oldv2 = oldhe2->vertex();
    std::list<CEdge*> & oldledges1 = oldv1->edges();
    std::list<CEdge*>::iterator it =  std::find(oldledges1.begin(), oldledges1.end(), e);
    if (it != oldledges1.end()) {
        // std::cout << "erase!" << "\n";
        oldledges1.erase(it);
    }

    std::list<CEdge*> & oldledges2 = oldv2->edges();
    it =  std::find(oldledges2.begin(), oldledges2.end(), e);
    if (it != oldledges2.end()) {
        // std::cout << "erase!" << "\n";
        oldledges2.erase(it);
    }
### adding a new edge
    CVertex* pV = (v1->id() < v2->id()) ? v1 : v2;
    std::list<CEdge*> & ledges = pV->edges();
    ledges.push_back(olded);
