#include <fstream>

#include "mesh/Mesh.h"

#include "mesh/RefDataRegular.h"
#include "linalg/Vector.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace mesh{

    Mesh::Mesh()
    {}
    
    Mesh::Mesh(const Mesh& other)
      : Nodes(other.Nodes), Cells(other.Cells)
    {
      copy(other);
    }

    Mesh& Mesh::operator=(const Mesh& other)
    {
      Nodes = other.Nodes;
      Cells = other.Cells;
      
      copy(other);

      return *this;
    }

    void Mesh::copy(const Mesh& other)
    {
      std::vector<Cell>::iterator it_new_cells;
      std::vector<Cell>::const_iterator it_cells;
      
      // Update pointers to local nodes
      for(it_cells = other.Cells.begin(), it_new_cells = Cells.begin();
	  it_cells != other.Cells.end(); ++it_cells, ++it_new_cells)
	{
	  for(int k=0; k<3; ++k)
	    it_new_cells->LocNode[k] = &(Nodes[it_cells->LocNode[k]->index]);	    
	}

      // \todo Find faster implementation. Edge list already created by the Refine routine.
      CreateEdgeList();
    }
    
    size_t Mesh::NrCells() const { return Cells.size(); }
    size_t Mesh::NrNodes() const { return Nodes.size(); }    

    void Mesh::WriteVtk(const std::string& filename) const
    {
      Vector x(NrNodes());
      WriteVtk(filename, x);
    }
    
    void Mesh::WriteVtk(const std::string& filename, const Vector& x) const
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
      size_t i;
      for (i=0, it = x.begin(); i<nr_nodes; ++it, ++i)
	ofs << *it << std::endl;
    }

    const std::vector<Cell>& Mesh::GetCellList() const
    {
      return Cells;
    }
    
    void Mesh::CreateEdgeList()
    {
      std::vector<Cell>::iterator cell;
      for(cell = Cells.begin(); cell != Cells.end(); ++cell)
	{
	  for(int i=0; i<3; ++i)
	    {
	      Node *Node0 = cell->LocNode[i], *Node1 = cell->LocNode[(i+1)%3];
	      Edge NewEdge(*cell, *Node0, *Node1);
	      
	      std::pair<std::set<Edge>::iterator, bool> it_pair = Edges.insert(NewEdge);
	      const Edge& CurEdge = *(it_pair.first);
	      
	      // If edge already exists update the neighbor
	      if (! it_pair.second)
		CurEdge.SetNeighbor(*cell);		

	      cell->LocEdge[i] = &CurEdge;
	    }
	}
    }

    Mesh& Mesh::RefineUniform()
    {
      std::vector<bool> cell_marker(NrCells(), true);            
      return Refine(cell_marker);
    }

    Mesh& Mesh::Refine(const std::vector<bool>& cell_marker)
    {
      Mesh* new_mesh = new Mesh();
      
      RefData ref_data = RefDataRegular();
      std::vector<Cell>::iterator it_cells;
      int j;
      for(it_cells = Cells.begin(), j=0;
	  it_cells != Cells.end(); ++j, ++it_cells)	
	it_cells->index = j;
      
      std::vector<Node*> new_edge_nodes(Edges.size(), NULL);  
      std::vector<const Edge*> new_edge_edges(2*Edges.size(), NULL);
      
      // Initialize edge index list (used to check if two cells share one new node)
      std::vector<size_t> global_edge_index(3*NrCells(), -1);

      size_t k;
      std::set<Edge>::const_iterator it_edge;
      
      for(k=0, it_edge = Edges.begin();
	  it_edge != Edges.end(); ++it_edge, ++k)
	{
	  for(int m=0; m<(it_edge->Type() == INTERFACE_EDGE ? 2 : 1); ++m)
	    {
	      Cell& cell = it_edge->GetNeighbor(m);
	      
	      size_t loc_edge_index = cell.EdgeIndex(*it_edge);
	      size_t cell_index = cell.Index();
	      
	      global_edge_index.at(3*cell_index + loc_edge_index) = k;
	    }
	}
      
      // Copy old nodes
      // Set index of node
      std::vector<Node>::iterator it_node;
      for(k=0, it_node = Nodes.begin(); it_node != Nodes.end(); ++it_node, ++k)
	it_node->index = k;
      
      // Perform copy operation
      for(size_t k=0; k<Nodes.size(); ++k)
	new_mesh->Nodes.push_back(Nodes.at(k));
      //new_mesh->Nodes.assign(Nodes.begin(), Nodes.end());

      // Initialize index counters
      size_t node_ctr = new_mesh->Nodes.size();
      size_t cell_ctr = 0;
	
      std::vector<bool>::const_iterator it_marker;
      std::vector<Cell>::const_iterator it_cell;
      for(j=0, it_marker = cell_marker.begin(), it_cell = Cells.begin();
	  it_cell != Cells.end(); ++it_cell, ++it_marker, ++j)
	{
	  if(!(*it_marker)) continue;
	  
	  // Regular refinement
	  int nr_nodes = ref_data.GetNrNodes();
	  int nr_edges = ref_data.GetNrEdges();

	  int edge_ctr = 0;
	  
	  std::vector<Node*> new_nodes(nr_nodes, NULL);
	  std::vector<const Edge*> new_edges(nr_edges, NULL);
	  
	  // Set old nodes
	  for(int k=0; k<3; ++k)	    
	    new_nodes.at(k) = &(new_mesh->Nodes.at(it_cell->LocNode[k]->Index()));	    
	  
	  // Create new nodes and edges
	  for(int k=0; k<3; ++k)
	    {
	      if(ref_data.GetEdgeVertex(k) == -1) continue;

	      size_t edge_index = global_edge_index.at(3*j + k);
	      if(new_edge_nodes.at(edge_index) == NULL)
		{
		  // Node does not exist. Create a new one.

		  // Compute vertex coordinates
		  double x = 0.5*(new_nodes.at(k)->getX() + new_nodes.at((k+1)%3)->getX());
		  double y = 0.5*(new_nodes.at(k)->getY() + new_nodes.at((k+1)%3)->getY());

		  // Create node object and insert into new mesh
		  new_mesh->Nodes.push_back(Node(node_ctr++, x, y));
		  
		  new_nodes.at(ref_data.GetEdgeVertex(k)) = &(new_mesh->Nodes.back());
		  new_edge_nodes.at(global_edge_index.at(3*j + k))
		    = new_nodes.at(ref_data.GetEdgeVertex(k));

		  // Create new edges
		  for(int m=0; m<2; ++m)
		    {
		      set_insert_res ins_res;
		      if(m==0)
			ins_res = new_mesh->Edges.insert(Edge(*new_nodes.at(k),
							      *new_nodes.at(ref_data.GetEdgeVertex(k))));
		      else
			ins_res = new_mesh->Edges.insert(Edge(*new_nodes.at(ref_data.GetEdgeVertex(k)),
							      *new_nodes.at((k+1)%3)));
		      
		      if(!ins_res.second)
			{
			  std::cerr << "Insertion of new edge " << m << " failed. Edge already exists.\n";			  
			  return *new_mesh;
			}
		      
		      new_edges.at(edge_ctr) = &(*ins_res.first);
		      new_edge_edges.at(2*edge_index+m) = new_edges.at(edge_ctr);

		      edge_ctr++;
		    }
		}
	      else
		{
		  // Node and edge already exist. Reuse the existing ones.
		  new_nodes.at(ref_data.GetEdgeVertex(k)) = new_edge_nodes.at(edge_index);

		  new_edges.at(edge_ctr++) = new_edge_edges.at(2*edge_index+1);
		  new_edges.at(edge_ctr++) = new_edge_edges.at(2*edge_index);		  
		}
	    }

	  // Create new interior edges
	  int nr_outer_edges = 2*ref_data.GetNrRefinedEdges();
	  if(nr_outer_edges != edge_ctr)
	    {
	      std::cerr << "Something went wrong. Number of outer edges is not correct.\n";
	      return *new_mesh;
	    }
	  
	  for(int k=nr_outer_edges; k<ref_data.GetNrEdges(); ++k)
	    {
	      const int* edge_vertices = ref_data.GetEdge(k);
	      set_insert_res ins_res = new_mesh->Edges.insert(Edge(*(new_nodes.at(edge_vertices[0])),
								   *(new_nodes.at(edge_vertices[1]))));
	      
	      new_edges.at(k) = &(*ins_res.first);
	    }
	  
	  // Create new cells
	  int nr_cells = ref_data.GetNrCells();
	  std::vector<Cell*> new_cells(nr_cells, NULL);
	  
	  for(int k=0; k<nr_cells; ++k)
	    {
	      const int* NewLocalNodes = ref_data.GetCell(k);

	      new_mesh->Cells.push_back(Cell(*(new_nodes.at(NewLocalNodes[0])),
					     *(new_nodes.at(NewLocalNodes[1])),
					     *(new_nodes.at(NewLocalNodes[2]))));
	      
	      new_cells.at(k) = &(new_mesh->Cells.back());
	      new_cells.at(k)->SetIndex(cell_ctr++);
	      
	      // Set pointers between edges and cells
	      const int* cell_edges = ref_data.GetCellEdges(k);
	      for(int m=0; m<3; ++m)
		{
		  new_edges.at(cell_edges[m])->SetNeighbor(*(new_cells.at(k)));
		  new_cells.at(k)->LocEdge[m] = new_edges.at(cell_edges[m]);
		}
	    }
	}
      
      return *new_mesh;
    }
    
    Mesh::~Mesh()
    {
      /// \todo Delete all Cells, Nodes and Edges.
    }
    
  };
};
