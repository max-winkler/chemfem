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
      size_t nr_dof = NrNodes();
      
      ofs << "# vtk DataFile Version 3.0" << std::endl;
      ofs << "2D scalar" <<std::endl;
      ofs << "ASCII" <<std::endl;
      ofs << "DATASET UNSTRUCTURED_GRID" <<std::endl;

      ofs << "POINTS " << nr_dof << " double" <<std::endl;
      for (std::vector<Node>::const_iterator it = Nodes.begin(); it != Nodes.end(); ++it)
	ofs << it->getX() << " " << it->getY() << " 0" << std::endl;
      
      ofs << "CELLS " << nr_cells << " " << 4*nr_cells << std::endl;
      for (std::vector<Cell>::const_iterator it = Cells.begin(); it != Cells.end(); ++it)	
	ofs << "3 " <<  it->LocNode[0]->Index << " " << it->LocNode[1]->Index
	    << " " << it->LocNode[2]->Index << std::endl;
      
      ofs << "CELL_TYPES " << nr_cells << std::endl;
      for (std::vector<Cell>::const_iterator it = Cells.begin(); it != Cells.end(); ++it)	
	ofs << 5 << std::endl;
      ofs << std::endl;

      // print data
      ofs <<  "POINT_DATA " << nr_dof << std::endl;
      ofs <<  "SCALARS nodes float" << std::endl;
      ofs <<  "LOOKUP_TABLE default"  << std::endl;

      // TODO: Works only for P1 elements
      for (Vector::const_iterator it = x.begin(); it != x.end(); ++it)
	ofs << *it << std::endl;
    }
    
  };
};
