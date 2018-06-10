#include <fstream>

#include "mesh/Mesh.h"

#include "linalg/Vector.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace mesh{

    size_t Mesh::NrCells() const { return Cells.size(); }
    size_t Mesh::NrNodes() const { return Nodes.size(); }    

    void Mesh::WriteVtk(const std::string& filename, const Vector& x)
    {
      std::ofstream ofs(filename);

      size_t nr_cells = NrCells();
      size_t nr_nodes = NrNodes();
      
      ofs << "# vtk DataFile Version 3.0" << std::endl;
      ofs << "2D scalar" <<std::endl;
      ofs << "ASCII" <<std::endl;
      ofs << "DATASET UNSTRUCTURED_GRID" <<std::endl;

      ofs << "POINTS " << nr_nodes << " double" <<std::endl;
      for (std::vector<Node>::const_iterator it = Nodes.begin(); it != Nodes.end(); ++it)
	ofs << it->getX() << " " << it->getY() << " 0" << std::endl;
      
      ofs << "CELLS " << nr_cells << " " << 4*nr_cells << std::endl;
      for (std::vector<Cell>::const_iterator it = Cells.begin(); it != Cells.end(); ++it)	
	ofs << "3 " <<  it->LocNode[0]->Index() << " " << it->LocNode[1]->Index()
	    << " " << it->LocNode[2]->Index() << std::endl;
      
      ofs << "CELL_TYPES " << nr_cells << std::endl;
      for (std::vector<Cell>::const_iterator it = Cells.begin(); it != Cells.end(); ++it)	
	ofs << 5 << std::endl;
      ofs << std::endl;

      // print data
      ofs <<  "POINT_DATA " << nr_nodes << std::endl;
      ofs <<  "SCALARS nodes float" << std::endl;
      ofs <<  "LOOKUP_TABLE default"  << std::endl;

      Vector::const_iterator it;
      int i;
      for (i=0, it = x.begin(); i<nr_nodes; ++it, ++i)
	ofs << *it << std::endl;
    }

    void Mesh::CreateEdgeList()
    {
      std::vector<Cell>::iterator cell;
      for(cell = Cells.begin(); cell != Cells.end(); ++cell)
	{
	  for(int i=0; i<3; ++i)
	    {
	      Node *Node0 = cell->LocNode[i], *Node1 = cell->LocNode[(i+1)%3];
	      Edge *NewEdge = new Edge(*cell, *Node0, *Node1);
	      
	      std::pair<std::set<Edge>::iterator, bool> it_pair = Edges.insert(*NewEdge);

	      // If edge already exists update the neighbor
	      if (! it_pair.second)
		{
		  const Edge& CurEdge = *(it_pair.first);
		  CurEdge.SetNeighbor(*cell);
		}
	    }
	}
    }

    Mesh::~Mesh()
    {
      // TODO: Delete all Cells, Nodes and Edges.
    }
    
  };
};
