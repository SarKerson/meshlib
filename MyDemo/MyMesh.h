 #ifndef _MY_MESH_
#define _MY_MESH_


#include "../MeshLib/core/Mesh/Vertex.h"
#include "../MeshLib/core/Mesh/Edge.h"
#include "../MeshLib/core/Mesh/Face.h"
#include "../MeshLib/core/Mesh/HalfEdge.h"
#include "../MeshLib/core/Mesh/BaseMesh.h"

#include "../MeshLib/core/Mesh/boundary.h"
#include "../MeshLib/core/Mesh/iterators.h"
#include "../MeshLib/core/Parser/parser.h"
#include <eigen3/Eigen/Dense>

#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

#define MESH_DEBUG 0

namespace MeshLib
{
	class CMyVertex;
	class CMyEdge;
	class CMyFace;
	class CMyHalfEdge;


	class CMyVertex : public CVertex
	{
	public:
		CMyVertex() : m_rgb(1, 1, 1) {};
		~CMyVertex() {};

		void _from_string();
		void _to_string();

		CPoint & rgb() { return m_rgb; };
	protected:
		CPoint m_rgb;
	};

	inline void CMyVertex::_from_string()
	{
		CParser parser(m_string);
		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken * token = *iter;
			if (token->m_key == "uv") //CPoint2
			{
				token->m_value >> m_uv;
			}
			if (token->m_key == "rgb") // CPoint
			{
				token->m_value >> m_rgb;
			}
		}
	}

	inline void CMyVertex::_to_string()
	{
		CParser parser(m_string);
		parser._removeToken("uv");

		parser._toString(m_string);
		std::stringstream iss;

		iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

		if (m_string.length() > 0)
		{
			m_string += " ";
		}
		m_string += iss.str();
	}

	class CMyEdge : public CEdge
	{
	public:
		CMyEdge() :m_sharp(false) {};
		~CMyEdge() {};

		void _from_string();
		void _to_string();

		bool & sharp() { return m_sharp; };
	protected:
		bool m_sharp;
	};

	inline void CMyEdge::_from_string()
	{
		CParser parser(m_string);
		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken * token = *iter;
			if (token->m_key == "sharp") // bool
			{
				m_sharp = true;
			}
		}
	}

	inline void CMyEdge::_to_string()
	{
		CParser parser(m_string);
		parser._removeToken("sharp");

		parser._toString(m_string);
		std::stringstream iss;

		if (m_sharp)
			iss << "sharp";

		if (m_string.length() > 0)
		{
			m_string += " ";
		}
		m_string += iss.str();
	}

	class CMyFace : public CFace
	{
	public:

		CPoint & normal() { return m_normal; };
	protected:
		CPoint m_normal;
	};

	class CMyHalfEdge : public CHalfEdge
	{
	};
    


	template<typename V, typename E, typename F, typename H>
	class MyMesh : public CBaseMesh<V, E, F, H>
	{
	public:
		
        // typedef V V;
        // typedef E E;
        // typedef F F;
        // typedef H H;

		typedef CBoundary<V, E, F, H>					CBoundary;
		typedef CLoop<V, E, F, H>						CLoop;

		typedef MeshVertexIterator<V, E, F, H>			MeshVertexIterator;
		typedef MeshEdgeIterator<V, E, F, H>			MeshEdgeIterator;
		typedef MeshFaceIterator<V, E, F, H>			MeshFaceIterator;
		typedef MeshHalfEdgeIterator<V, E, F, H>		MeshHalfEdgeIterator;

		typedef VertexVertexIterator<V, E, F, H>		VertexVertexIterator;
		typedef VertexEdgeIterator<V, E, F, H>			VertexEdgeIterator;
		typedef VertexFaceIterator<V, E, F, H>			VertexFaceIterator;
		typedef VertexInHalfedgeIterator<V, E, F, H>	VertexInHalfedgeIterator;
		typedef VertexOutHalfedgeIterator<V, E, F, H>	VertexOutHalfedgeIterator;

		typedef FaceVertexIterator<V, E, F, H>			FaceVertexIterator;
		typedef FaceEdgeIterator<V, E, F, H>			FaceEdgeIterator;
		typedef FaceHalfedgeIterator<V, E, F, H>		FaceHalfedgeIterator;

		void output_mesh_info();
		void test_iterator();

	};

	typedef MyMesh<CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> CMyMesh;

	template<typename V, typename E, typename F, typename H>
	void MeshLib::MyMesh<V, E, F, H>::output_mesh_info()
	{
		int nv = this->numVertices();
		int ne = this->numEdges();
		int nf = this->numFaces();

		std::cout << "#V=" << nv << "  ";
		std::cout << "#E=" << ne << "  ";
		std::cout << "#F=" << nf << "  ";

		int euler_char = nv - ne + nf;
		std::cout << "Euler's characteristic=" << euler_char << "  ";

		CBoundary boundary(this);
		std::vector<CLoop*> & loops = boundary.loops();
		int nb = loops.size();

		int genus = (2 - (euler_char + nb)) / 2;
		std::cout << "genus=" << genus << std::endl;
	}


	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::test_iterator()
	{
		for (MeshVertexIterator viter(this); !viter.end(); ++viter)
		{
			V * pV = *viter;
			// you can do something to the vertex here
			// ...

			for (VertexVertexIterator vviter(pV); !vviter.end(); ++vviter)
			{
				V * pW = *vviter;
				// you can do something to the neighboring vertices with CCW
				// ...
			}

			for (VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
			{
				E * pE = *veiter;
				// you can do something to the neighboring edges with CCW
				// ...
			}

			for (VertexFaceIterator vfiter(pV); !vfiter.end(); ++vfiter)
			{
				F * pF = *vfiter;
				// you can do something to the neighboring faces with CCW
				// ...
			}

			for (VertexInHalfedgeIterator vhiter(this, pV); !vhiter.end(); ++vhiter)
			{
				H * pH = *vhiter;
				// you can do something to the incoming halfedges with CCW
				// ...
			}
		}

		for (MeshEdgeIterator eiter(this); !eiter.end(); ++eiter)
		{
			E * pE = *eiter;
			// you can do something to the edge here
			// ...
		}

		for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter)
		{
			F * pF = *fiter;
			// you can do something to the face here
			// ...
		}

		//there are some other iterators which you can find them in class MyMesh

		std::cout << "Iterators test OK.\n";
	}
}

inline double mod(MeshLib::CPoint v) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline double mod(Eigen::Vector3d v) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


namespace gauss_bonnet {

	using namespace MeshLib;

    /*
    My codes to test THE GAUSS-BONNET THEOROM:
    */
    double checkG_B(CMyMesh *pMesh)
    {
        double GB = 0;
        for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
        {
            CMyVertex *v = *viter;
            if(v->boundary()) {     //pi
                std::vector<CPoint> vEdge;
                double totalAngle = 0;
                for(CMyMesh::VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
                {
                    CMyVertex *pV = *vviter;
                    CPoint p = pV->point();
                    CPoint e = p - v->point();   //e1
                    vEdge.push_back(e);
                }
                for(size_t i = 1; i < vEdge.size(); ++i) {
                    CPoint pa = vEdge.at(i - 1),
                           pb = vEdge.at(i);
                    double cos_ = pa * pb / mod(pa) / mod(pb);
                    totalAngle += acos(cos_);
                }
                GB += (M_PI - totalAngle);
            }
            if(!v->boundary()) {     //2 * pi
                std::vector<CPoint> vEdge;
                double totalAngle = 0;
                for(CMyMesh::VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
                {
                    CMyVertex *pV = *vviter;
                    CPoint p = pV->point();
                    CPoint e = p - v->point();   //e1
                    vEdge.push_back(e);
                }
                for(size_t i = 1; i < vEdge.size(); ++i) {
                    CPoint pa = vEdge.at(i - 1),
                           pb = vEdge.at(i);
                    double cos_ = pa * pb / mod(pa) / mod(pb);
                    totalAngle += acos(cos_);
                }
                //below is the important codes without witch the result will be wrong
                //the angle between the fist vector and the last vector
                CPoint p_first = vEdge.at(0), p_last = vEdge.at(vEdge.size() - 1);
                double cos_ = p_first * p_last / mod(p_first) / mod(p_last);
                totalAngle += acos(cos_);
                GB += (2 * M_PI - totalAngle);
            }
            
        }
        return GB;
    }

}


namespace convexHull {
	using namespace MeshLib;

    typedef struct myEdge 
    {
        CVertex* _vertex;
        CVertex* _source;
        myEdge* next;
        myEdge* prev;


        myEdge() 
        {
            _vertex = NULL;
            _source = NULL;
            next = NULL;
            prev = NULL;
        }

        CVertex*& vertex() { return _vertex; }
        CVertex*& source() { return _source; }
        myEdge*& he_next() { return next; }
        myEdge*& he_prev() { return prev; }

    } myEdge;

    void my_read_m(CMyMesh &mesh, const char *input ) {
        std::fstream is( input, std::fstream::in );

        if( is.fail() )
        {
            fprintf(stderr,"Error in opening file %s\n", input );
            return;
        }

        char buffer[MAX_LINE];
        int id;

        while( is.getline(buffer, MAX_LINE )  )
        {       
        
            std::string line( buffer );
            line = strutil::trim( line );

            strutil::Tokenizer stokenizer( line, " \r\n" );

            stokenizer.nextToken();
            std::string token = stokenizer.getToken();
        
            if( token == "Vertex"  ) 
            {
                stokenizer.nextToken();
                token = stokenizer.getToken();
                id = strutil::parseString<int>(token);

                CPoint p;
                for( int i = 0 ; i < 3; i ++ )
                {
                    stokenizer.nextToken();
                    token = stokenizer.getToken();
                    p[i] = strutil::parseString<float>(token);
                }
            
                CVertex* v  = mesh.createVertex( id );
                v->point() = p;
                v->id()    = id;

                if( ! stokenizer.nextToken("\t\r\n") ) continue;
                token = stokenizer.getToken();

                int sp = (int) token.find("{");
                int ep = (int) token.find("}");

                if( sp >= 0 && ep >= 0 )
                {
                    v->string() = token.substr( sp+1, ep-sp-1 );
                }
                continue;
            }
            

            if( token == "Face" ) 
            {

                stokenizer.nextToken();
                token = stokenizer.getToken();
                id = strutil::parseString<int>(token);
        
                std::vector<CMyVertex*> v;

                while( stokenizer.nextToken() )
                {
                    token = stokenizer.getToken();
                    if( strutil::startsWith( token, "{" ) ) break;
                    int vid = strutil::parseString<int>(token);
                    v.push_back( mesh.idVertex( vid ) );
                }

                CFace *f = mesh.createFace( v, id );

                if( ! stokenizer.nextToken("\t\r\n") ) continue;
                token = stokenizer.getToken();

                //stokenizer.reset();
                token = line;
                int sp = (int) token.find("{");
                int ep = (int) token.find("}");

                if( sp >= 0 && ep >= 0 )
                {
                    f->string() = token.substr( sp+1, ep-sp-1 );
                }
                continue;
            }
        }
    }

    std::vector<myEdge*> lboundary;

    void printBoundaryInfo() {
        std::cout << "convexhull: ";
        myEdge* edge = lboundary[0], *t = edge->he_next();
        std::cout << edge->vertex()->point() << " ->\t";
        while (t->vertex() != edge->source()) {
            std::cout << t->vertex()->point() << " ->\t";
            t = t->he_next();
        }
        std::cout << t->vertex()->point() << "\n";
        std::cout << "total : " << lboundary.size() << "\n";
    }

    void printFaceInfo(CFace* f)
    {
        CHalfEdge *he = f->halfedge();
        std::cout << "info of face: \n\tp1: " << he->vertex()->point() << "\n"             //i
        << "\tp2: " << he->he_next()->vertex()->point() << "\n"                            //j
        << "\tp3: " << he->he_prev()->vertex()->point() << "\n";                           //k
    }

    double orientation(CPoint pi, CPoint pj, /*CPoint *pk,*/ CPoint p)   //i->j->p
    {
        Eigen::Matrix2d m;
        m << pj[0] - pi[0], p[0] - pi[0],
             pj[1] - pi[1], p[1] - pi[1];
        return m.determinant();
    }

    double orientation(myEdge* he, CVertex* v) 
    {
        double result = orientation(he->source()->point(),
                           he->vertex()->point(),
                           v->point());
        // std::cout << "from " << he->source()->point() << " to " << he->vertex()->point() << " with " << v->point() << "\t:" << result << "\n";
        return result;
    }

    void connectBoundary(CVertex* v) {   //all is new pointer!!!!
        myEdge *next = NULL, *pre = NULL;
        for (std::vector<myEdge*>::iterator it = lboundary.begin(); it!= lboundary.end(); ++it) {
            if ((*it)->he_next() == NULL) {
                pre = *it;
            }
            if ((*it)->he_prev() == NULL) {
                next = *it;
            }
        }
        assert ( pre != NULL && next != NULL);

        myEdge *_nhe1 = new myEdge(),
               *_nhe2 = new myEdge();

        _nhe1->vertex() = v;
        _nhe1->source() = pre->vertex();
        _nhe2->vertex() = next->source();
        _nhe2->source() = v;
        //linking to each other
        pre->he_next() = _nhe1;
        _nhe1->he_prev() = pre;
        _nhe1->he_next() = _nhe2;
        _nhe2->he_prev() = _nhe1;
        _nhe2->he_next() = next;
        next->he_prev() = _nhe2;

        lboundary.push_back(_nhe1);
        lboundary.push_back(_nhe2);
    }

    void makeConvexHull(CMyMesh &mesh) {
        CFace *face = mesh.idFace(1);
        myEdge    *h0 = new myEdge(),
                  *h1 = new myEdge(),
                  *h2 = new myEdge();
        h0->vertex() = face->halfedge()->vertex();
        h1->vertex() = face->halfedge()->he_next()->vertex();
        h2->vertex() = face->halfedge()->he_prev()->vertex();
        lboundary.push_back(h0);
        lboundary.push_back(h1);
        lboundary.push_back(h2);
        //linking to each other
        for( size_t i = 0; i < lboundary.size(); i ++ )
        {
            lboundary[i]->he_next() = lboundary[(i+1)%3];
            lboundary[i]->he_prev() = lboundary[(i+2)%3];
            lboundary[i]->source() = lboundary[i]->he_prev()->vertex();
        }

        for (CMyMesh::MeshVertexIterator viter(&mesh); !viter.end(); ++viter) {
            //picked one point

            CVertex* _v_ = *viter;
            if (_v_ != h0->vertex()
                && _v_ != h1->vertex()
                && _v_ != h2->vertex()) {
                //for each edge in boundary
                bool flag = false;
                for (std::vector<myEdge*>::iterator it = lboundary.begin(); it!= lboundary.end();) {
                    myEdge *_he_ = *it;
                    double result = orientation(_he_, _v_);
                    if (result < 0) {         //!!!!!!!!!!!!!!!
                        flag = true;
                        myEdge* &he = _he_;
                        it = lboundary.erase(it);
                        if (he->he_prev())
                            he->he_prev()->he_next() = NULL;
                        if (he->he_next())
                            he->he_next()->he_prev() = NULL;  //so as to find boundary
                        delete he;
                    } else {
                        ++it;
                    }
                }
                if (flag) {
                    connectBoundary(_v_);
                }
            }
        }
    }

}

namespace delaunayTriangulation
{
	using namespace MeshLib;
    void printFaceInfo(CFace* f)
    {
        CHalfEdge *he = f->halfedge();
        std::cout << "info of face: \n\tp1: " << he->vertex()->point() << "\n"             //i
        << "\tp2: " << he->he_next()->vertex()->point() << "\n"                            //j
        << "\tp3: " << he->he_prev()->vertex()->point() << "\n";                           //k
    }

    double inCircle(CPoint pi, CPoint pj, CPoint pk, CPoint p)
    {
        Eigen::Matrix4d m;
        m << pi[0], pi[1], pi[0] * pi[0] + pi[1] * pi[1], 1,
             pj[0], pj[1], pj[0] * pj[0] + pj[1] * pj[1], 1,
             pk[0], pk[1], pk[0] * pk[0] + pk[1] * pk[1], 1,
             p[0],  p[1],  p[0] *  p[0] +  p[1] *  p[1], 1;
        return m.determinant();
    }


    double triArea(CPoint p0, CPoint p1, CPoint p2)
    {
        Eigen::Matrix3d m;
        m << p0[0], p0[1], 1,
             p1[0], p1[1], 1,
             p2[0], p2[1], 1;
        return m.determinant() * 0.5;
    }


    CPoint triBaryCentric(CPoint pi, CPoint pj, CPoint pk, CPoint p)
    {
        CPoint result;
        double area = triArea(pi, pj, pk);
        result[0] = triArea(p, pi, pj) / area;
        result[1] = triArea(p, pj, pk) / area;
        result[2] = triArea(p, pk, pi) / area;
        return result;
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


    CFace* LocatePoint(CMyMesh &mesh, CPoint p)
    {
        CFace *face = mesh.idFace(1);  //initial face
        while (true) {
            CHalfEdge *he = face->halfedge();
            CPoint barycentric = triBaryCentric(he->vertex()->point(),              //i
                                                he->he_next()->vertex()->point(),   //j
                                                he->he_prev()->vertex()->point(),   //k
                                                p);
            if (barycentric[0] < 0) {   //ij
                he = he->he_next()->he_sym();
                if (he != NULL) {
                    face = he->face();
                    continue;
                }
            }
            he = face->halfedge();        //important code!!!'casuse the halfedge might change
            if (barycentric[1] < 0) {    //jk
                he = he->he_prev()->he_sym();
                if (he != NULL) {
                    face = he->face();
                    continue;
                }
            }
            he = face->halfedge();
            if (barycentric[2] < 0) {
                he = he->he_sym();
                if (he != NULL) {
                    face = he->face();
                    continue;
                }
            }
            if (barycentric[0] > 0 && barycentric[1] > 0 && barycentric[2] > 0) {
                return face;
            }
            else {
                return NULL;
            }
        }
    }


    void faceSplit(CMyMesh &mesh, CFace* f, CPoint p)
    {
        int ivid = mesh.vertices().size() + 1;
        CVertex *v = mesh.createVertex(ivid);
        v->point() = p;
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

    }


    bool legalizeEdge(CMyMesh &mesh, CVertex* v, CEdge* e)
    {
        if (e->boundary()) {
            return false;
        }

        CHalfEdge* h = e->halfedge(0);
        CVertex *v0, *v1;
        if (h->he_next()->vertex() == v) {
            assert(h != NULL);
            h = e->halfedge(1);
        }
        v0 = h->vertex();
        v1 = h->he_prev()->vertex();

        CVertex* v2 = h->he_next()->vertex();
        double result = inCircle(v0->point(), v1->point(), v->point(), v2->point());
        if (result > 0) {
            edgeSwap((CMyEdge*)e);
            legalizeEdge(mesh, v, mesh.createEdge((CMyVertex*)v1, (CMyVertex*)v2));
            legalizeEdge(mesh, v, mesh.createEdge((CMyVertex*)v0, (CMyVertex*)v2));
            return true;
        }
        else {
            return false;
        }
    }

    void insertVertex(CMyMesh& mesh, CPoint p)
    {
        CFace* face = LocatePoint(mesh, p);
        if (face != NULL) {
            // printFaceInfo(face);
            CHalfEdge *h0 = face->halfedge(),
                      *h1 = h0->he_next(),
                      *h2 = h0->he_prev();
            faceSplit(mesh, face, p);
            assert(h0->he_next()->vertex() == h1->he_next()->vertex() 
                   && h0->he_next()->vertex() ==h2->he_next()->vertex());
            CVertex* v = h0->he_next()->vertex();
            legalizeEdge(mesh, v, h0->edge());
            legalizeEdge(mesh, v, h1->edge());
            legalizeEdge(mesh, v, h2->edge());                //problem!!!!!!!!!!!!
        } else {
            std::cout << p << "---out---\n";
        }
    }

    double miniFaceAngle(CFace * face) {

    }

    void legalizeFace(CFace * face) {

    }

}

namespace harmonicMap {

    using namespace MeshLib;
    double Energy = 0;
    double e = 0.01;

/**
 * function of the boudary vertexs, [x0, y0, z0],[x1, y1, z1],...,[xn, yn, zn]
 */
    std::map<CMyVertex*, Eigen::Vector3d> g;    //all the boundary

    std::map<CMyVertex*, Eigen::Vector3d> f;    //all the result 



    float k(CMyVertex * v, CMyVertex * w, CMyMesh & mesh)
    {
        CEdge * edge =  mesh.vertexEdge(v, w);
        Eigen::Vector3d v_v; v_v << (v->point())[0], (v->point())[1], (v->point())[2];
        Eigen::Vector3d v_w; v_w << (w->point())[0], (w->point())[1], (w->point())[2];
        assert(edge != NULL);
        if (edge->boundary()) {
            CHalfEdge * he = (edge->halfedge(0))->he_next();
            assert(he != NULL);
            CVertex * vk = he->vertex();
            Eigen::Vector3d v_k; v_k << (vk->point())[0], (vk->point())[1], (vk->point())[2];
            double cot_k = (v_v - v_k).dot(v_w - v_k) / (mod((v_v - v_k).cross(v_w - v_k)));
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
            Eigen::Vector3d v_k; v_k << (vk->point())[0], (vk->point())[1], (vk->point())[2];
            Eigen::Vector3d v_l; v_l << (vl->point())[0], (vl->point())[1], (vl->point())[2];
            /* calculate the cot */
#if MESH_DEBUG
            std::cout << "vv - vk" << v_v - v_k << "\n"
                      << "vw - vk" << v_w - v_k << "\n"
                      << "vv - vk dot vw - vk: " << (v_v - v_k).dot(v_w - v_k) << "\n";
#endif
            double cot_k = (v_v - v_k).dot(v_w - v_k) / (mod((v_v - v_k).cross(v_w - v_k)));
            double cot_l = (v_v - v_l).dot(v_w - v_l) / (mod((v_v - v_l).cross(v_w - v_l)));   
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

    void initDisk(CMyMesh & mesh)
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

        for(int i = 0; i < _vlist_.size(); ++i) {
            double theta = vlength[i] / total_length * 2.0 * M_PI;
            Eigen::Vector3d g_v;
            g_v << cos(theta), sin(theta), 0.0; 
            g.insert(std::pair<CMyVertex*, Eigen::Vector3d>(_vlist_[i], g_v));
            f.insert(std::pair<CMyVertex*, Eigen::Vector3d>(_vlist_[i], g_v));
            _vlist_[i]->point() = CPoint(g_v[0], g_v[1], g_v[2]);
#if MESH_DEBUG
            std::cout << g_v << "\n";
#endif
        }
    }

    void topologicalDisk(CMyMesh & mesh)
    {
        initDisk(mesh);
        int num_boundary = num_of_boundary(mesh);
        int num_inner = num_of_inner(mesh);
#if MESH_DEBUG
        std::cout << "num_boundary: " << num_boundary << "\n"
                  << "num_inner: " << num_inner << "\n";
#endif

        assert(num_boundary == g.size());        //all the functions of boundary vertex

        double E = 0;
        if (num_inner != 0) {                     //inner vertex

            std::map<CMyVertex*, int> index;

            int i = 0;
            for (CMyMesh::MeshVertexIterator viter(&mesh); !viter.end(); ++viter)
            {
                CMyVertex *v = *viter;
                if (!v->boundary()) {            //all inner vertex
                    index.insert(std::make_pair(v, i++));
                }
            }

            for (std::map<CMyVertex*, int>::iterator it = index.begin(); it != index.end(); ++ it) {
                CMyVertex* v = it->first;
                Eigen::Vector3d value; value << 0.0, 0.0, 0.0;
                f.insert(std::make_pair(v, value));
            }                                    //init f, now for each vertex in the Mesh, f(v) = [0 , 0, 0]
            

            i = 0;
            do {                                        

                Energy = E;
                E = 0;

                Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_inner, num_inner);
                Eigen::VectorXd bx(num_inner),
                                by(num_inner);

                assert(index.size() == num_inner);
                for (std::map<CMyVertex*, int>::iterator it = index.begin(); it != index.end(); ++ it) {
                    CMyVertex* v = it->first;
                    for(CMyMesh::VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
                    {
                        CMyVertex *w = *vviter;
                        double k_v_w = k(v, w, mesh);
                        if (!w->boundary()) { //both inner
                            A(index[v], index[w]) = k_v_w;
                            A(index[w], index[v]) = k_v_w;
                            A(index[v], index[v]) -= k_v_w;
                            A(index[w], index[w]) -= k_v_w;  
                        }
                        else {
                            A(index[v], index[v]) -= k_v_w;
                            bx(index[v]) -= k_v_w * g[w][0];
                            by(index[v]) -= k_v_w * g[w][1];
                        }
                        E += k_v_w * (f[v] - f[w]).dot(f[v] - f[w]);
                    }
                }
                std::cout << "E: " << E << "\n";
                std::cout << "Energy: " << Energy << "\n";
                std::cout << "diff: " << fabs(E - Energy) << "\n";
                // Ax = bx, Ay = by
                // get x, y
                Eigen::VectorXd x = A.colPivHouseholderQr().solve(bx);
                Eigen::VectorXd y = A.colPivHouseholderQr().solve(by);

                for (std::map<CMyVertex*, int>::iterator it = index.begin(); it != index.end(); ++ it) {
                    CMyVertex* v = it->first;
                    int indx = it->second;
                    Eigen::Vector3d value; value << x[indx], y[indx], 0.0;
                    f[v] = value;
                    v->point() = CPoint(value[0], value[1], value[2]);
                }
                ++i;
            } while(fabs(Energy - E) < e || i < 10);

        }// if

        // for (std::map<CMyVertex*, Eigen::Vector3d>::iterator it = g.begin(); it != g.end(); ++ it) {
        //     f[it->first] = it->second;
        // }

        for (std::map<CMyVertex*, Eigen::Vector3d>::iterator it = f.begin(); it != f.end(); ++ it) {
            Eigen::Vector3d v = it->second;
            it->first->point() = CPoint(v[0], v[1], v[2]);   //change all the vertex
#if MESH_DEBUG
            std::cout << "result: " << it->first->point() << "\n"; 
#endif
        }
    }



}

#endif // !_MY_MESH_
