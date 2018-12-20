#include "mesh/RefData.h"

namespace chemfem{
  namespace mesh{

    RefData::RefData(int NrNewCells, int NrNewNodes, int NrNewEdges)
      : NrNewCells(NrNewCells), NrNewNodes(NrNewNodes), NrNewEdges(NrNewEdges)
    {
      // NewNodeCoords = new double[2*NrNewNodes];
      // NewCells = new int[3*NrNewCells];
      // NewEdges = new int[2*NrNewEdges];
    }

    RefData::~RefData()
    {
      delete[] NewNodeCoords;
      delete[] NewCells;
      delete[] NewEdges;
    }
    
    int RefData::GetNrNodes()
    { return NrNewNodes; }

    int RefData::GetNrEdges()
    { return NrNewEdges; }
    
    int RefData::GetNrCells()
    { return NrNewCells; }

    const double* RefData::GetNodeCoords(int i)
    { return NewNodeCoords[i]; }

    const int* RefData::GetEdge(int i)
    { return NewEdges[i]; }

    const int* RefData::GetCell(int i)
    { return NewCells[i]; }
    
    const int RefData::GetEdgeVertex(int i)
    { return NewEdgeVertices[i]; }
  }
}
  
