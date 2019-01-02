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
      std::cout << "Copy constructor called\n";
    }

    Mesh& Mesh::operator=(const Mesh& other)
    {
      copy(other);
      std::cout << "Copy assignment called\n";      
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
      std::vector<size_t> global_edge_index(3*NrCells());

      size_t k;
      std::set<Edge>::const_iterator it_edge;
      
      for(k=0, it_edge = Edges.begin();
	  it_edge != Edges.end(); ++it_edge, ++k)
	{
	  for(int m=0; m<(it_edge->Type() == INTERFACE_EDGE ? 2 : 1); ++m)
	    {
	      Cell& cell = it_edge->GetNeighbor(m);
	      global_edge_index[3*cell.Index() + cell.EdgeIndex(*it_edge)] = k;
	    }
	}

      // Copy old nodes
      // Set index of node
      std::vector<Node>::iterator it_node;
      for(k=0, it_node = Nodes.begin(); it_node != Nodes.end(); ++it_node, ++k)
	it_node->index = k;
      
      // Perform copy operation
      new_mesh->Nodes.assign(Nodes.begin(), Nodes.end());

      size_t node_ctr = new_mesh->Nodes.size();
	
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
	  
	  Node** new_nodes = new Node*[nr_nodes];
	  const Edge** new_edges = new const Edge*[nr_edges];

	  for(int k=0; k<nr_nodes; ++k)
	    new_nodes[k] = NULL;
	  for(int k=0; k<nr_edges; ++k)
	    new_edges[k] = NULL;
	  
	  // Set old nodes
	  for(int k=0; k<3; ++k)	    
	    new_nodes[k] = &(new_mesh->Nodes[it_cell->LocNode[k]->Index()]);	    
	  
	  // Create new nodes and edges
	  for(int k=0; k<3; ++k)
	    {
	      if(ref_data.GetEdgeVertex(k) == -1) continue;

	      size_t edge_index = global_edge_index[3*j + k];
	      if(new_edge_nodes[edge_index] == NULL)
		{
		  // Node does not exist. Create a new one.

		  // Compute vertex coordinates
		  double x = 0.5*(new_nodes[k]->getX() + new_nodes[(k+1)%3]->getX());
		  double y = 0.5*(new_nodes[k]->getY() + new_nodes[(k+1)%3]->getY());

		  // Create node object and insert into new mesh
		  new_mesh->Nodes.push_back(Node(node_ctr++, x, y));
		  
		  new_nodes[ref_data.GetEdgeVertex(k)] = &(new_mesh->Nodes.back());
		  new_edge_nodes[global_edge_index[3*j + k]] = new_nodes[ref_data.GetEdgeVertex(k)];

		  // Create new edges
		  for(int m=0; m<2; ++m)
		    {
		      set_insert_res ins_res;
		      if(m==0)
			ins_res = new_mesh->Edges.insert(Edge(*new_nodes[k],
							     *new_nodes[ref_data.GetEdgeVertex(k)]));
		      else
			ins_res = new_mesh->Edges.insert(Edge(*new_nodes[ref_data.GetEdgeVertex(k)],
							     *new_nodes[(k+1)%3]));
		      
		      if(!ins_res.second)
			{
			  std::cerr << "Insertion of new edge " << m << " failed. Edge already exists.\n";
			  std::cerr << "Trying to add edges between " << *new_nodes[k] << " and "
				    << *new_nodes[ref_data.GetEdgeVertex(k)] << std::endl;
			  			  
			  for(std::set<Edge>::const_iterator it_edges = new_mesh->Edges.begin();
			      it_edges != new_mesh->Edges.end(); ++it_edges)
			    {
			      std::cout << "Edge: " << *it_edges << std::endl;
			    }
			  
			  return *new_mesh;
			}
		      
		      new_edges[edge_ctr] = &(*ins_res.first);
		      new_edge_edges[2*edge_index+m] = new_edges[edge_ctr];

		      edge_ctr++;
		    }
		}
	      else
		{
		  // Node and edge already exist. Reuse the existing ones.
		  new_nodes[ref_data.GetEdgeVertex(k)] = new_edge_nodes[edge_index];

		  new_edges[edge_ctr++] = new_edge_edges[2*edge_index+1];
		  new_edges[edge_ctr++] = new_edge_edges[2*edge_index];		  
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
	      set_insert_res ins_res = new_mesh->Edges.insert(Edge(*(new_nodes[edge_vertices[0]]),
							       *(new_nodes[edge_vertices[1]])));
	      
	      new_edges[k] = &(*ins_res.first);
	    }
	  
	  // Create new cells
	  int nr_cells = ref_data.GetNrCells();
	  Cell** new_cells = new Cell*[nr_cells];
	  
	  for(int k=0; k<nr_cells; ++k)
	    {
	      const int* NewLocalNodes = ref_data.GetCell(k);

	      new_mesh->Cells.push_back(Cell(*(new_nodes[NewLocalNodes[0]]),
					    *(new_nodes[NewLocalNodes[1]]),
					    *(new_nodes[NewLocalNodes[2]])));
	      
	      new_cells[k] = &(new_mesh->Cells.back());

	      // Set pointers between edges and cells
	      const int* cell_edges = ref_data.GetCellEdges(k);
	      for(int m=0; m<3; ++m)
		{
		  new_edges[cell_edges[m]]->SetNeighbor(*(new_cells[k]));
		  new_cells[k]->LocEdge[m] = new_edges[cell_edges[m]];
		}
	    }	 	  
	}

      return *new_mesh;
    }
    
    Mesh::~Mesh()
    {
      std::cout << "Destructor called\n";
      /// \todo Delete all Cells, Nodes and Edges.
    }
    
  };
};
