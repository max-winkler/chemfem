#include "mesh/RefData.h"

namespace chemfem{
  namespace mesh{

    RefData::RefData(int NrNewCells, int NrNewNodes, int NrNewEdges)
      : NrNewCells(NrNewCells), NrNewNodes(NrNewNodes), NrNewEdges(NrNewEdges)
    { }
    
    
    int RefData::GetNrNodes()
    { return NrNewNodes; }

    int RefData::GetNrEdges()
    { return NrNewEdges; }
    
    int RefData::GetNrCells()
    { return NrNewCells; }

    const double* RefData::GetNodeCoords(int i)
    { return &(NewNodeCoords[3*i]); }

    const int* RefData::GetEdge(int i)
    { return &(NewEdges[2*i]); }

    const int* RefData::GetCell(int i)
    { return &(NewCells[3*i]); }

  }
}

