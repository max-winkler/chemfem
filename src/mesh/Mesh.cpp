#include <fstream>
#include <deque>
#include <iomanip>
#include <map>
#include <algorithm>

#include "mesh/Mesh.h"

#include "mesh/RefDataRegular.h"
#include "mesh/RefDataBisection0.h"
#include "mesh/RefDataBisection1.h"
#include "mesh/RefDataBisection2.h"

#include "linalg/Vector.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace mesh{

    Mesh::Mesh()
    {}
    
    Mesh::Mesh(const Mesh& other)
      : Nodes(other.Nodes), Cells(other.Cells), Edges(other.Edges)
    {
      std::cout << "Mesh copied\n";
      copy(other);
    }

    Mesh& Mesh::operator=(const Mesh& other)
    {
      std::cout << "Mesh copied by =\n";
      Nodes = other.Nodes;
      Cells = other.Cells;
      Edges = other.Edges;
      
      copy(other);
      
      return *this;
    }

    void Mesh::copy(const Mesh& other)
    {
      // \todo Find faster implementation. Edge list already created by the Refine routine.
      // CreateEdgeList();
    }
    
    std::ostream& operator<<(std::ostream& os, const Mesh& mesh)
    {
      os << "Vertices:\n";
      
      std::vector<Node>::const_iterator it_nodes;
      size_t k = 0;
      for(it_nodes = mesh.Nodes.begin(), k=0; it_nodes != mesh.Nodes.end(); ++it_nodes, ++k)
        os << " " << std::setw(3) << k << " : " << (*it_nodes) << std::endl;

      os << "\nCells:\n";
      std::vector<Cell>::const_iterator it_cells;
      for(it_cells = mesh.Cells.begin(), k=0; it_cells != mesh.Cells.end(); ++it_cells, ++k)
        {
	os << " " << std::setw(3) << k << " : " << (*it_cells) << std::endl;
        }

      os << "\nEdges:\n";
      std::vector<Edge>::const_iterator it_edges;
      for(it_edges = mesh.Edges.begin(), k=0; it_edges != mesh.Edges.end(); ++it_edges, ++k)
        {
	os << " " << std::setw(3) << k << " : " << (*it_edges) << std::endl;
        }
      return os;
    }
    
    size_t Mesh::NrCells() const { return Cells.size(); }
    size_t Mesh::NrNodes() const { return Nodes.size(); }
    size_t Mesh::NrEdges() const { return Edges.size(); }    

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
	ofs << "3 " <<  it->LocNode[0] << " " << it->LocNode[1]
	    << " " << it->LocNode[2] << std::endl;
      
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
      std::cerr << "WARNING: The routine CreateEdgeList is very slow and should be avoided\n";
      
      // Remove old edges
      Edges.clear();
      
      std::vector<Cell>::iterator cell;
      size_t cell_ind;
      
      // Remove all edge information
      for(cell = Cells.begin(), cell_ind=0; cell != Cells.end(); ++cell, ++cell_ind)
	{
	  cell->SetIndex(cell_ind);
	  for(int k=0; k<3; ++k)
	    cell->LocEdge[k] = -1;
	}

      for(cell = Cells.begin(); cell != Cells.end(); ++cell)
	{
	  for(int i=0; i<3; ++i)
	    {
	      size_t Node0 = cell->LocNode[i], Node1 = cell->LocNode[(i+1)%3];
	      Edge new_edge(Node0, Node1, BOUNDARY_EDGE);

	      // Check if edge already exists
	      std::vector<Edge>::iterator it = std::find(Edges.begin(), Edges.end(), new_edge);
	      if(it != Edges.end())
	        {
		// Edge exists, update the neighbor
		(*it).SetNeighbor(cell->Index());
		cell->LocEdge[i] = std::distance(Edges.begin(), it);
	        }
	      else
	        {
		// Edge does not exist, add to edge list
		new_edge.SetNeighbor(cell->Index());
		Edges.push_back(new_edge);
		cell->LocEdge[i] = Edges.size()-1;
	        }
	    }
	}
    }

    Mesh& Mesh::RefineUniform()
    {
      std::vector<bool> cell_marker(NrCells(), true);            
      return Refine(cell_marker);
    }

    struct EdgeRefInfo
    {
      size_t new_edge[2];
      size_t new_inner_vertex;
      /**
      EdgeRefInfo() {}
      EdgeRefInfo(const EdgeRefInfo& info) : new_inner_vertex(info.new_inner_vertex)
      {
        std::cout << "EDGE REF COPIED\n";
        new_edge[0] = info.new_edge[0];
        new_edge[1] = info.new_edge[1];
      }
      */
    };
    
    Mesh& Mesh::Refine(const std::vector<bool>& cell_marker)
    {
      // Copy the mesh
      Mesh* new_mesh = new Mesh(*this);

      std::cout << "New mesh:\n" << *new_mesh << std::endl;

      // Refinement descriptors
      RefData ref_data[]
        = {RefDataBisection0(), RefDataBisection1(), RefDataBisection2()};
      RefDataBisection0 tmp;
      
      std::deque<size_t> marked_cell_idx;

      std::map<size_t, EdgeRefInfo> edge_ref_infos;
      size_t k = 0;

      // Initialize stack of marked cells
      std::vector<bool>::const_iterator it_marker;
      for(it_marker = cell_marker.begin(), k=0; it_marker != cell_marker.end(); ++it_marker, ++k)
        if(*it_marker) marked_cell_idx.push_back(k);
      
      while(!marked_cell_idx.empty())
        {
	size_t idx = marked_cell_idx.front();
	marked_cell_idx.pop_front();
	
	std::cout << "I am about to refine element " << idx << std::endl;

	Cell cur_cell = new_mesh->Cells[idx];

	// determine largest index
	size_t max_idx = 0;
	if(cur_cell.LocNode[1] > cur_cell.LocNode[0]) max_idx = 1;
	if(cur_cell.LocNode[2] > cur_cell.LocNode[max_idx]) max_idx = 2;
	
	// Print edge ref infos
	std::cout << "Edge ref infos available for ";
	for(std::map<size_t, EdgeRefInfo>::const_iterator it = edge_ref_infos.begin();
	    it != edge_ref_infos.end(); ++it)
	  {
	    std::cout << it->first << " ";
	  }
	std::cout << std::endl;
	
	// Store edge index
	size_t ref_edge_idx = cur_cell.LocEdge[(max_idx+1)%3];
	bool ref_edge_refined = (edge_ref_infos.find(ref_edge_idx) != edge_ref_infos.end());
	Edge ref_edge = new_mesh->Edges[ref_edge_idx];
	
	EdgeRefInfo edge_ref_info;
	if(ref_edge_refined)
	  edge_ref_info = edge_ref_infos[ref_edge_idx];
	
	std::cout << "I will refine edge " << ref_edge_idx << std::endl;
	
	// apply refinement to the cell
	RefData& cur_ref_data = ref_data[max_idx];
	
	// Create new vertices
	size_t new_nodes[4];

	for(size_t i=0; i<3; ++i)
	  new_nodes[i] = cur_cell.LocNode[i];
          	  
	if(ref_edge_refined)
	  {
	    // Use vertex already created by neighbour
	    new_nodes[3] = edge_ref_info.new_inner_vertex;
	    std::cout << "EDGE REFINEMENT INFO FOUND! We use here node " << new_nodes[3] << std::endl;
	  }
	else
	  {
	    // Calculate vertex coordinates
	    double s = cur_ref_data.NewNodeCoords[0][0];
	    double t = cur_ref_data.NewNodeCoords[0][1];

	    double new_x = (1-s-t)*Nodes[new_nodes[0]].getX()
	      + s*Nodes[new_nodes[1]].getX()
	      + t*Nodes[new_nodes[2]].getX();
	    double new_y = (1-s-t)*Nodes[new_nodes[0]].getY()
	      + s*Nodes[new_nodes[1]].getY()
	      + t*Nodes[new_nodes[2]].getY();

	    // Add vertex to vertex list
	    new_mesh->Nodes.push_back(Node(new_mesh->NrNodes(), new_x, new_y));
	    new_nodes[3] = new_mesh->Nodes.size()-1;
	  }

	// Create new edges
	int nr_new_edges = cur_ref_data.GetNrEdges();
	size_t* new_edges = new size_t[nr_new_edges];

	for(int i=0; i<nr_new_edges; ++i)
	  new_edges[i] = -1;         
	
	for(int i=0; i<3; ++i)
	  {
	    // Copy old edges 
	    if(cur_ref_data.OldEdgeNewEdgeLen[i] == 1)
	      {
	        new_edges[cur_ref_data.OldEdgeNewEdge[i][0]] = cur_cell.LocEdge[i];

	        // Unset edge-cell relation
	        if(edge_ref_infos.find(cur_cell.LocEdge[i]) == edge_ref_infos.end())
		{
		  new_mesh->Edges[cur_cell.LocEdge[i]].UnsetNeighbor(idx);
		}
	      }	   	      
	    else
	      {
	        if(ref_edge_refined)
		{
		  // New edges were created by neighbor already
		  for(int j=0; j<cur_ref_data.OldEdgeNewEdgeLen[i]; ++j)		    
		    new_edges[cur_ref_data.OldEdgeNewEdge[i][j]] = edge_ref_info.new_edge[1-j];
		}
	        else
		{		  	        
		  // Create new edges
		  for(int j=0; j<cur_ref_data.OldEdgeNewEdgeLen[i]; ++j)
		    {
		      int new_edge_idx = cur_ref_data.OldEdgeNewEdge[i][j];

		      Edge edge(new_nodes[cur_ref_data.NewEdges[new_edge_idx][0]],
			      new_nodes[cur_ref_data.NewEdges[new_edge_idx][1]],
			      EdgeType::BOUNDARY_EDGE);
		  
		      if(j == 0)
		        {
			new_mesh->Edges[cur_cell.LocEdge[i]] = edge;
			new_edges[new_edge_idx] = cur_cell.LocEdge[i];
		        }
		      else
		        {
			new_mesh->Edges.push_back(edge);		     
			new_edges[new_edge_idx] = new_mesh->Edges.size()-1;
		        }
		    }
		}
	      }
	  }

	// Create interior edges
	for(int i=0; i<cur_ref_data.InteriorEdgesLen; ++i)
	  {
	    Edge edge(new_nodes[cur_ref_data.NewEdges[cur_ref_data.InteriorEdges[i]][0]],
		    new_nodes[cur_ref_data.NewEdges[cur_ref_data.InteriorEdges[i]][1]],
		    EdgeType::BOUNDARY_EDGE);
	    
	    new_mesh->Edges.push_back(edge);
	    new_edges[cur_ref_data.InteriorEdges[i]] = new_mesh->Edges.size()-1;
	  }

	size_t* new_cells = new size_t[cur_ref_data.NrNewCells];
	
	// Create new cells
	for(int i=0; i<cur_ref_data.GetNrCells(); ++i)
	  {
	    // Get nodes of new cell
	    size_t loc_nodes[3];
	    for(int j=0; j<3; ++j)
	        loc_nodes[j] = new_nodes[cur_ref_data.NewCells[i][j]];

	    // Create new cell
	    Cell new_cell(loc_nodes[0], loc_nodes[1], loc_nodes[2]);

	    // Create cell-edge relations
	    for(int j=0; j<3; ++j)
	      new_cell.LocEdge[j] = new_edges[cur_ref_data.NewCellEdges[i][j]];
	    
	    // Insert cell into cell list
	    if(i==0)
	      {
	        // Overwrite parent cell	        
	        new_mesh->Cells[idx] = new_cell;
	        new_cells[i] = idx;
	      }
	    else
	      {
		// Create new cell at end of the list
		new_mesh->Cells.push_back(new_cell);
		new_cells[i] = new_mesh->Cells.size()-1;
	      }
	  }

	// Create edge-cell relations for current cell
	for(int i=0; i<cur_ref_data.NrNewCells; ++i)
	  {
	    for(int j=0; j<3; ++j)
	      {
	        int edge_idx = cur_ref_data.NewCellEdges[i][j];

	        // Do not set neighbor when edge was refined by neighbor
	        if(edge_ref_infos.find(new_edges[edge_idx]) != edge_ref_infos.end())
		continue;
	        
	        new_mesh->Edges[new_edges[edge_idx]].SetNeighbor(new_cells[i]);
	      }
	  }

	// Add child to marked cell list
	for(int i=0; i<3; ++i)
	  {
	    if(cur_ref_data.OldEdgeNewEdgeLen[i] > 1) continue;
	    
	    if(edge_ref_infos.find(cur_cell.LocEdge[i]) != edge_ref_infos.end())
	      {
	        // Get index of new cell
	        size_t cell_idx = new_cells[cur_ref_data.OldEdgeNewCell[i]];
	        std::cout << "I also have to refine the child cell " << cell_idx << std::endl;
	        marked_cell_idx.push_back(cell_idx);
	      }
	  }
	
	// Add neighbor to marked cell list
	if(ref_edge.Type() == INTERFACE_EDGE)
	  {
	    size_t neigh_cell = ref_edge.GetNeighbor(idx);
	    if(std::find(marked_cell_idx.begin(), marked_cell_idx.end(), neigh_cell) != marked_cell_idx.end())
	      continue;

	    std::cout << "I also have to refine the neighbor with index " << neigh_cell << std::endl;
	    marked_cell_idx.push_back(neigh_cell);

	    // Create refinement information for the neighbor	        		      
	    EdgeRefInfo ref_info;
	    ref_info.new_inner_vertex = new_nodes[cur_ref_data.NewEdgeVertices[(max_idx+1)%3]];
	    for(int k=0; k<2; ++k)		
	      ref_info.new_edge[k] = new_edges[cur_ref_data.OldEdgeNewEdge[(max_idx+1)%3][k]];

	    std::cout << "CREATED Refinement info for edge " << ref_edge_idx << std::endl;
	    std::cout << " inner node is " << ref_info.new_inner_vertex
		    << ", new edges are " << new_edges[cur_ref_data.OldEdgeNewEdge[(max_idx+1)%3][0]]
		    << " and " << new_edges[cur_ref_data.OldEdgeNewEdge[(max_idx+1)%3][1]] << std::endl;
	    edge_ref_infos[ref_edge_idx] = ref_info;
	  }

	// Delete edge ref info
	if(ref_edge_refined)
	  edge_ref_infos.erase(ref_edge_idx);
	
	std::cout << "mesh after refinement step:\n" << *new_mesh << std::endl;
        }     

      // not good but for testing only
      //new_mesh->CreateEdgeList();
      std::cout << "final mesh:\n" << *new_mesh << std::endl;
      
      return *new_mesh;
    }
    
    /*
      // RED-GREEN Refinement (not working yet)
    Mesh& Mesh::Refine(const std::vector<bool>& cell_marker)
    {
      Mesh* new_mesh = new Mesh();

      // Reserve memory for Nodes and Cells
      new_mesh->Cells.reserve(4*NrCells());
      new_mesh->Nodes.reserve(NrNodes()+NrEdges());      
      new_mesh->BdEdges.reserve(2*BdEdges.size());
      
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
      // \todo: The following method is probably faster
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
	  
	  // Create new edges
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

		      // If the edge is a boundary edge, at it to the BdEdges list
		      if(it_cell->LocEdge[k]->Type() == BOUNDARY_EDGE)
		        new_mesh->BdEdges.push_back(new_edges.at(edge_ctr));
		      
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
						         *(new_nodes.at(edge_vertices[1]))
						         ));
	      
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
	      
	      new_cells[k] = &(new_mesh->Cells.back());
	      new_cells[k]->SetIndex(cell_ctr++);

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
    */
    
    Mesh::~Mesh()
    {
      /// \todo Delete all Cells, Nodes and Edges.
    }

    bool Mesh::Check()
    {
      for(std::vector<Cell>::const_iterator it_cell = Cells.begin(); it_cell != Cells.end(); ++it_cell)
	{
	  const Cell& cell = *it_cell;

	  for(int k=0; k<3; ++k)
	    {
	      const size_t edge_idx = cell.LocEdge[k];

	      if(cell.EdgeIndex(edge_idx) != k)
	        std::cerr << "The edge index seems to be wrong.\n";
	      
	      if(!((Edges[edge_idx].Node0 == cell.LocNode[k]
		    && Edges[edge_idx].Node1 == cell.LocNode[(k+1)%3])
		   || (Edges[edge_idx].Node1 == cell.LocNode[k]
		       && Edges[edge_idx].Node0 == cell.LocNode[(k+1)%3])))
		{
		  std::cerr << "I found a cell (index " << it_cell->Index() << ") "
			    << "whose edge (index " << edge_idx
			    << ") does not point to the corresponding nodes.\n";
		  return false;
		}
	    }
	}
      
      return true;      
    }

    double Mesh::Determinant(size_t i) const
    {
      double x[3], y[3];
      for(int k=0; k<3; ++k)
	{
	  x[k] = Nodes[Cells[i].LocNode[k]].getX();
	  y[k] = Nodes[Cells[i].LocNode[k]].getY();
	}
      
      return (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]);
    }

    using chemfem::linalg::DenseMatrix;
    
    DenseMatrix Mesh::Jacobian(size_t i) const
    {
      double x[3], y[3];
      for(int k=0; k<3; ++k)
	{
	  x[k] = Nodes[Cells[i].LocNode[k]].getX();
	  y[k] = Nodes[Cells[i].LocNode[k]].getY();
	}
      
      DenseMatrix Jac(2,2);
      Jac.data[0] =  x[1] - x[0];
      Jac.data[1] =  x[2] - x[0];
      Jac.data[2] =  y[1] - y[0];
      Jac.data[3] =  y[2] - y[0];

      return Jac;
    }

  };
};

