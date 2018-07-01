#include "mesh/RefData.h"

namespace chemfem{
  namespace mesh{

    RefData::RefData(int NrNewCells, int NrNewNodes, int NrNewEdges)
      : NrNewCells(NrNewCells), NrNewNodes(NrNewNodes), NrNewEdges(NrNewEdges)
    {
      NewNodeCoords = new double[2*NrNewNodes];
      NewCells = new int[3*NrNewCells];
      NewEdges = new int[2*NrNewEdges];
    }
    
    
    int RefData::GetNrNodes()
    { return NrNewNodes; }

    int RefData::GetNrEdges()
    { return NrNewEdges; }
    
    int RefData::GetNrCells()
    { return NrNewCells; }

    double* RefData::GetNodeCoords(int i)
    { return &(NewNodeCoords[3*i]); }

    int* RefData::GetEdge(int i)
    { return &(NewEdges[2*i]); }

    int* RefData::GetCell(int i)
    { return &(NewCells[3*i]); }

  }
}
  
