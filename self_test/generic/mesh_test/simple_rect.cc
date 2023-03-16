#include <iostream>

#include "generic.h"
#include "meshes.h"

using namespace std;
using namespace oomph;

class QFaceElement : public virtual FaceElement,
                     // public virtual FaceGeometry<QElement<2, 2>>,
                     public virtual QElement<1, 2>
{
};

template<class ELEMENT>
class FaceMesh : public Mesh
{
public:
  /*FaceMesh(Mesh* const& bulk_mesh_pt, unsigned const& face_boundary_id) :
  Mesh()
  {
    if (face_boundary_id > bulk_mesh_pt->nboundary())
    {
      cout << "error" << endl;
    }
    else
    {
      Vector<unsigned> other_boundary_id;
      for (unsigned b = 0; b < bulk_mesh_pt->nboundary(); b++)
      {
        if (b != face_boundary_id)
        {
          other_boundary_id.push_back(b);
        }
      }
      unsigned boundary_count = 0;

      // Find the number of nodes on the boundary
      unsigned nbound_node = bulk_mesh_pt->nboundary_node(face_boundary_id);
      // Loop over the boundary nodes and add them to face mesh node pointer
      for (unsigned n = 0; n < nbound_node; n++)
      {
        BoundaryNode* nod_pt = dynamic_cast<BoundaryNode*>(
          bulk_mesh_pt->boundary_node_pt(face_boundary_id, n));
        bool isBoundaryNode = false;
        // Check if it is another boundary node
        for (unsigned i = 0; i < other_boundary_id.size(); i++)
        {
          if (nod_pt->is_on_boundary(other_boundary_id[i]))
          {
          }
        }
        if (isBoundaryNode)
        {
          add_boundary_node(3, nod_pt);
        }
        else
        {
          // otherwise, add as an interior node
          this->add_node_pt(nod_pt);
        }
      }

      // Find the number of elements next to the boundary
      unsigned nbound_element =
        bulk_mesh_pt->nboundary_element(face_boundary_id);

      // Allocate the store for the elements
      Element_pt.resize(nbound_element);

      // Loop over the elements adjacent to boundary b
      for (unsigned e = 0; e < nbound_element; e++)
      {
        Element_pt[e] = new ELEMENT;

        // Save bulk element and face index
        // Create the FaceElement
        Bulk_element_pt.push_back(
          bulk_mesh_pt->boundary_element_pt(face_boundary_id, e));
        Bulk_face_index.push_back(
          bulk_mesh_pt->face_index_at_boundary(face_boundary_id, e));
      }
    }
  }

  virtual ~FaceMesh()
  {
    // Free the elements
    // Loop over the elements in reverse order
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long i = Element_pt_range; i > 0; i--)
    {
      delete Element_pt[i - 1];
      Element_pt[i - 1] = 0;
    }

    // Assume the bulk mesh will delete the nodes we have added
    Node_pt.clear();
  }

private:
  Vector<GeneralisedElement*> Bulk_element_pt;
  Vector<int> Bulk_face_index;
  */
};

class FaceMesh<SimpleRectangularQuadMesh> : public LineMeshBase
{
};


int main()
{
  cout << "Mesh test" << endl;
  const unsigned Nx = 4;
  const unsigned Ny = 4;
  const double Lx = 1;
  const double Ly = 1;
  SimpleRectangularQuadMesh<QElement<2, 2>> mesh(Nx, Ny, Lx, Ly);

  cout << mesh.nboundary() << endl;
  const unsigned face_boundary_id = 2;
  FaceMesh<QFaceElement> face_mesh(&mesh, face_boundary_id);

  cout << face_mesh.nboundary() << endl;

  return 0;
}
