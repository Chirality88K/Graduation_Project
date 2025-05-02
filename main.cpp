#include "opennurbs.h"
#include "example_ud.h"
#include "read.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "FilletSurface.h"
#include "PiecesClothoid3d.h"
#include "Cone_Surface.h"
#include "Spiral.h"
#include "Cornu_spiral.h"
#include <algorithm>
#include "EulerBezier2d.hpp"
const double PI = acos(-1.0);


using namespace std;

//static bool write_cylinder(const wchar_t* filename, ON_TextLog& error_log);
//static bool write_circle(const wchar_t* filename, ON_TextLog& error_log);
void Add_Cylinder(ONX_Model* model, const wchar_t* name, double radius = 1, double height = 2, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Add_Circle(ONX_Model* model, const wchar_t* name, double radius = 1, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Add_Ball(ONX_Model* model, const wchar_t* name, double radius = 1, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Add_BezierCurve(ONX_Model* model, const wchar_t* name, const ON_BezierCurve* bc, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Get_Curvature_and_Torsion_of_NurbsCurve(const ON_NurbsCurve& onc, string filename);
void Cornu_test(ONX_Model* model);
void EulerBsplineSpiralInterpolation(ON_NurbsCurve& onc, int max_vtx_num = 50);
void EulerBsplineTest(ONX_Model* model);

ON_3dmObjectAttributes* Internal_CreateManagedAttributes(
  int layer_index,
  const wchar_t* name
)
{
  ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
  attributes->m_layer_index = layer_index;
  attributes->m_name = name;
  return attributes;
}

static bool write_points_example( const wchar_t* filename, ON_TextLog& error_log  )
{
  // example demonstrates how to write a singe points and point clouds
  ONX_Model model;
  INTERNAL_INITIALIZE_MODEL(model);

  // file settings (units, tolerances, views, ...)
  // OPTIONAL - change values from defaults
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::LengthUnitSystem::Meters;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.01;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%

  // layer table
  // define some layers
  model.AddDefaultLayer(nullptr, ON_Color::UnsetColor);
  const int point1_layer_index = model.AddLayer(L"my layer",ON_Color::Black);
  const int pointcloud_layer_index = model.AddLayer(L"red points",ON_Color::SaturatedRed);
  const int point2_layer_index = model.AddLayer(L"one blue point",ON_Color::SaturatedBlue);

  // we'll put 2 red and one blue point in a group
  ON_Group group;
  group.SetName(L"group of points");
  group.SetIndex(0);
  model.AddModelComponent(group, true);


  // single point at (1,4,5) on default layer
  ON_Point* point1 = new ON_Point(ON_3dPoint( 1.0, 4.0, 5.0 ));
  point1->AttachUserData( new CExampleWriteUserData("write_points_example()-point1") );
  model.AddManagedModelGeometryComponent( 
    point1,
    Internal_CreateManagedAttributes(point1_layer_index,L"first point")
  );

  // point "cloud" with 3 points on red point cloud layer
  ON_PointCloud* pointcloud = new ON_PointCloud();
  pointcloud->AppendPoint(ON_3dPoint( 1.0, 6.0, 5.0 ));
  pointcloud->AppendPoint(ON_3dPoint( 1.5, 4.5, 6.0 ));
  pointcloud->AppendPoint(ON_3dPoint( 2.0, 5.0, 7.0 ));
  pointcloud->AttachUserData( new CExampleWriteUserData("write_points_example()-pointcloud") );
  ON_3dmObjectAttributes* pointcloud_attributes = Internal_CreateManagedAttributes(pointcloud_layer_index, L"3 points");
  pointcloud_attributes->AddToGroup(group.Index());
  model.AddManagedModelGeometryComponent( 
    pointcloud,
    pointcloud_attributes
  );

  // single point at (3,2,4) on red point layer and in group with the pointcloud
  ON_Point* point2 = new ON_Point(ON_3dPoint( 3.0, 2.0, 4.0  ));
  ON_3dmObjectAttributes* point2_attributes = Internal_CreateManagedAttributes(point2_layer_index, L"last point");
  point2_attributes->AddToGroup(group.Index());
  point2->AttachUserData( new CExampleWriteUserData("write_points_example()-point2") );
  model.AddManagedModelGeometryComponent( point2, point2_attributes);

  return Internal_WriteExampleModel(model, filename, error_log);
}

static bool write_curves_example( const wchar_t* filename, ON_TextLog& error_log )
{
  // example demonstrates how to write a NURBS curve, line, and circle
  ONX_Model model;
  INTERNAL_INITIALIZE_MODEL(model);
  
  // file settings (units, tolerances, views, ...)
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::LengthUnitSystem::Inches;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%

  // add some layers
  model.AddDefaultLayer(nullptr, ON_Color::UnsetColor);
  const int line_layer_index = model.AddLayer(L"line layer",ON_Color::Black);
  const int wiggle_layer_index = model.AddLayer(L"green NURBS wiggle",ON_Color::SaturatedGreen);
  const int circles_layer_index = model.AddLayer(L"blue circles",ON_Color::SaturatedBlue);
  
  {
    // add a line
    ON_Object* managed_line = new ON_LineCurve( ON_Line( ON_3dPoint(1.0,2.0,-1.5), ON_3dPoint(5.0,3.0,2.0) ) );
    model.AddManagedModelGeometryComponent(
      managed_line,
      Internal_CreateManagedAttributes(line_layer_index, L"straight line curve")
      );
  }

  {
    // add a wiggly cubic curve
    ON_NurbsCurve* wiggle = new ON_NurbsCurve(
      3, // dimension
      false, // true if rational
      4,     // order = degree+1
      6      // number of control vertices
      );
    int i;
    for ( i = 0; i < wiggle->CVCount(); i++ ) {
      ON_3dPoint pt( 2*i, -i, (i-3)*(i-3) ); // pt = some 3d point
      wiggle->SetCV( i, pt );
    }

    // ON_NurbsCurve's have order+cv_count-2 knots.
    wiggle->SetKnot(0, 0.0);
    wiggle->SetKnot(1, 0.0);
    wiggle->SetKnot(2, 0.0);
    wiggle->SetKnot(3, 1.5);
    wiggle->SetKnot(4, 2.3);
    wiggle->SetKnot(5, 4.0);
    wiggle->SetKnot(6, 4.0);
    wiggle->SetKnot(7, 4.0);

    model.AddManagedModelGeometryComponent(
      wiggle,
      Internal_CreateManagedAttributes(wiggle_layer_index, L"wiggly cubic curve")
      );
  }

  {
    // add two circles
    ON_ArcCurve* circle1 = new ON_ArcCurve( ON_Circle( ON_3dPoint(1.0,2.0,-1.5), 3.0 ) );
    model.AddManagedModelGeometryComponent(
      circle1,
      Internal_CreateManagedAttributes(circles_layer_index, L"radius 3 circle")
      );

    ON_ArcCurve* circle2 = new ON_ArcCurve( ON_Circle( ON_3dPoint(1.0,2.0,-1.5), 5.0 ) );
    model.AddManagedModelGeometryComponent(
      circle2,
      Internal_CreateManagedAttributes(circles_layer_index, L"radius 5 circle")
      );
  }

  return Internal_WriteExampleModel(model, filename, error_log);
}


static bool write_surfaces_example( const wchar_t* filename, ON_TextLog& error_log )
{
  // example demonstrates how to write a NURBS surface
  ONX_Model model;
  INTERNAL_INITIALIZE_MODEL(model);

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////

  // The code between the comment bands has nothing to do with I/O.
  // It is simply an easy way to get a NURBS surface to write.
  const int bIsRational = false;
  const int dim = 3;
  const int u_degree = 2;
  const int v_degree = 3;
  const int u_cv_count = 3;
  const int v_cv_count = 5;

  // The knot vectors do NOT have the 2 superfluous knots
  // at the start and end of the knot vector.  If you are
  // coming from a system that has the 2 superfluous knots,
  // just ignore them when writing a 3dm file.
  double u_knot[ u_cv_count + u_degree - 1 ];
  double v_knot[ v_cv_count + v_degree - 1 ];

  // make up a quadratic knot vector with no interior knots
  u_knot[0] = u_knot[1] = 0.0;
  u_knot[2] = u_knot[3] = 1.0;

  // make up a cubic knot vector with one simple interior knot
  v_knot[0] = v_knot[1] = v_knot[2] = 0.0;
  v_knot[3] = 1.5;
  v_knot[4] = v_knot[5] = v_knot[6] = 2.0;

  // Rational control points can be in either homogeneous
  // or euclidean form. Non-rational control points do not
  // need to specify a weight.  
  ON_3dPoint CV[u_cv_count][v_cv_count];

  int i, j;
  for ( i = 0; i < u_cv_count; i++ ) {
    for ( j = 0; j < v_cv_count; j++ ) {
      CV[i][j].x = i;
      CV[i][j].y = j;
      CV[i][j].z = i-j;
    }
  }

  // write a line on the default layer
  ON_NurbsSurface nurbs_surface( dim, bIsRational, 
                        u_degree+1, v_degree+1,
                        u_cv_count, v_cv_count );

  for ( i = 0; i < nurbs_surface.KnotCount(0); i++ )
    nurbs_surface.SetKnot( 0, i, u_knot[i] );

  for ( j = 0; j < nurbs_surface.KnotCount(1); j++ )
    nurbs_surface.SetKnot( 1, j, v_knot[j] );

  for ( i = 0; i < nurbs_surface.CVCount(0); i++ ) {
    for ( j = 0; j < nurbs_surface.CVCount(1); j++ ) {
      nurbs_surface.SetCV( i, j, CV[i][j] );
    }
  }

  model.AddModelGeometryComponent(&nurbs_surface, nullptr);
  //   model.AddDefaultLayer(L"NURBS surface", ON_Color::UnsetColor);

  return Internal_WriteExampleModel(model, filename, error_log);
}


static bool write_mesh_example( const wchar_t* filename, ON_TextLog& error_log )
{
  // example demonstrates how to create and write a mesh
  ONX_Model model;
  INTERNAL_INITIALIZE_MODEL(model);

  model.AddDefaultLayer(L"mesh", ON_Color::Black);

  // create a mesh to write
  // The mesh is a pyramid with 4 triangular sides and a quadranglar 
  // base.  The mesh has 5 vertices and 5 faces.  
  // The side faces share normals at their common vertices.  The
  // quadrangular base has normals different from the side normal.
  // Coincident vertices that have distinct normals must be
  // duplicated in the vertex list.
  //
  // The apex will be at (1,1.5,4) with normal (0,0,1).
  // The base corners will be at (0,0,0), (0,2,0), (2,3,0), (0,3,0).


  bool bHasVertexNormals = true; // we will specify vertex normals
  bool bHasTexCoords = false;    // we will not specify texture coordinates
  const int vertex_count = 5+4;  // 4 duplicates for different base normals
  const int face_count = 5; // 4 triangle sides and a quad base
  ON_Mesh mesh( face_count, vertex_count, bHasVertexNormals, bHasTexCoords);

  // The SetVertex(), SetNormal(), SetTCoord() and SetFace() functions
  // return true if successful and false if input is illegal.  It is
  // a good idea to inspect this returned value.

  // vertex #0: apex location and normal
  mesh.SetVertex( 0, ON_3dPoint(1.0,  1.5,  5.0) );
  mesh.SetVertexNormal( 0, ON_3dVector(0.0,  0.0,  1.0) );

  // vertex #1: SW corner vertex for sides
  mesh.SetVertex( 1, ON_3dPoint(0.0,  0.0,  0.0) );
  mesh.SetVertexNormal( 1, ON_3dVector(-1.0, -1.0,  0.0) ); // set normal will unitize if needed

  // vertex #2: SE corner vertex for sides
  mesh.SetVertex( 2, ON_3dPoint(2.0,  0.0,  0.0) );
  mesh.SetVertexNormal( 2, ON_3dVector(+1.0, -1.0,  0.0) );

  // vertex #3: NE corner vertex for sides
  mesh.SetVertex( 3, ON_3dPoint(2.0,  3.0,  0.0) );
  mesh.SetVertexNormal( 3, ON_3dVector(+1.0, +1.0,  0.0) );

  // vertex #4: NW corner vertex for sides
  mesh.SetVertex( 4, ON_3dPoint(0.0,  3.0,  0.0) );
  mesh.SetVertexNormal( 4, ON_3dVector(-1.0, +1.0,  0.0) );

  // vertex #5: SW corner vertex for base
  mesh.SetVertex( 5, ON_3dPoint(0.0,  0.0,  0.0) ); // == location of v1
  mesh.SetVertexNormal( 5, ON_3dVector(0.0,  0.0, -1.0) );

  // vertex #6: SE corner vertex for base
  mesh.SetVertex( 6, ON_3dPoint(2.0,  0.0,  0.0) ); // == location of v2
  mesh.SetVertexNormal( 6, ON_3dVector(0.0,  0.0, -1.0) );

  // vertex #7: SW corner vertex for base
  mesh.SetVertex( 7, ON_3dPoint(2.0,  3.0,  0.0) ); // == location of v3
  mesh.SetVertexNormal( 7, ON_3dVector(0.0,  0.0, -1.0) );

  // vertex #8: SW corner vertex for base
  mesh.SetVertex( 8, ON_3dPoint(0.0,  3.0,  0.0) ); // == location of v4
  mesh.SetVertexNormal( 8, ON_3dVector(0.0,  0.0, -1.0) );

  // faces have vertices ordered counter-clockwise

  // South side triangle
  mesh.SetTriangle( 0,   1, 2, 0 );

  // East side triangle
  mesh.SetTriangle( 1,   2, 3, 0 );

  // North side triangle
  mesh.SetTriangle( 2,   3, 4, 0 );

  // West side triangle
  mesh.SetTriangle( 3,   4, 1, 0 );

  // last face is quadrangular base
  mesh.SetQuad( 4,   5, 8, 7, 6 );

  if ( !mesh.HasVertexNormals() )
    mesh.ComputeVertexNormals();

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////

  // Avoid copying the mesh - useful technique for large objects
  model.AddModelGeometryComponentForExperts(false, &mesh, false, nullptr, true);

  return Internal_WriteExampleModel(model, filename, error_log);
}

static void make_trimming_curves( ON_Brep& brep, 
                                  const ON_2dPoint& A2, // start point in parameter space
                                  const ON_2dPoint& B2, // end point in parameter space
                                  const ON_3dPoint& A3, // start point in parameter space
                                  const ON_3dPoint& B3  // end point in parameter space
                                  )
{
  ON_LineCurve* p2dCurve = new ON_LineCurve( A2, B2 );
  ON_LineCurve* p3dCurve = new ON_LineCurve( A3, B3 );

  // it is not necessary for the domains of the 2d and 3d curves
  // to match, but it makes it easier to understand the brep
  ON_Interval domain = p3dCurve->Domain();
  p2dCurve->SetDomain( domain.Min(), domain.Max() );

  brep.m_C2.Append(p2dCurve);

  brep.m_C3.Append(p3dCurve);
}


static bool write_trimmed_surface_example( const wchar_t* filename, ON_TextLog& error_log )
{
  // write a trimmed surface
  ONX_Model model;
  INTERNAL_INITIALIZE_MODEL(model);

  model.AddDefaultLayer(L"trimmed surface", ON_Color::Black);

  // trimmed surfaces are written as a CRhinoBrep that has
  // a single surface and a single CRhinoBrepFace.
  //
  // Trimming loops are simple closed curves and are oriented
  // so that the active portion of the trimmed surface's
  // domain lies to the left of the trimming curves.

  ON_Brep brep;
  ON_2dPoint q;

  // Create a 10x10 plane surface at z=3 with domain [0,1]x[0,1]
  ON_PlaneSurface* pSurface = new ON_PlaneSurface( ON_Plane( ON_3dPoint( 0, 0,3), 
                                                             ON_3dPoint(10,10,3), 
                                                             ON_3dPoint(10, 0,3) ) );
  pSurface->SetDomain(0,0.0,10.0);
  pSurface->SetDomain(1,0.0,10.0);

  // ~ON_Brep() will delete this surface
  const int si = brep.m_S.Count(); // index of surface
  brep.m_S.Append(pSurface);

  // create simple trimming triangle
  ON_2dPoint A2(1.0, 2.0); // parameter space locations of 2d trim corners
  ON_2dPoint B2(9.0, 1.5);
  ON_2dPoint C2(7.0, 8.0);

  ON_3dPoint A3 = pSurface->PointAt(A2.x,A2.y);
  ON_3dPoint B3 = pSurface->PointAt(B2.x,B2.y);
  ON_3dPoint C3 = pSurface->PointAt(C2.x,C2.y);

  make_trimming_curves( brep, A2, B2, A3, B3 ); // creates 2d and 3d curve
  make_trimming_curves( brep, B2, C2, B3, C3 );
  make_trimming_curves( brep, C2, A2, C3, A3 );

  // there are vertices at the 3 corners
  brep.NewVertex( pSurface->PointAt( A2.x, A2.y ) );
  brep.NewVertex( pSurface->PointAt( B2.x, B2.y ) );
  brep.NewVertex( pSurface->PointAt( C2.x, C2.y ) );

  // the vertices are exact since we have lines on a plane
  brep.m_V[0].m_tolerance = 0.0;
  brep.m_V[1].m_tolerance = 0.0;
  brep.m_V[2].m_tolerance = 0.0;

  // there are 3 edges along the sides of the triangle
  brep.NewEdge( brep.m_V[0], brep.m_V[1], 0 ); // start vertex, end vertex, 3d curve index
  brep.NewEdge( brep.m_V[1], brep.m_V[2], 1 ); // start vertex, end vertex, 3d curve index
  brep.NewEdge( brep.m_V[2], brep.m_V[0], 2 ); // start vertex, end vertex, 3d curve index

  // the edges are exact since we have lines on a plane
  brep.m_E[0].m_tolerance = 0.0;
  brep.m_E[1].m_tolerance = 0.0;
  brep.m_E[2].m_tolerance = 0.0;

  // there is 1 face
  ON_BrepFace& face = brep.NewFace( si );

  // outer boundary trimming loops
  ON_BrepLoop& loop = brep.NewLoop( ON_BrepLoop::outer, face );

  // geometrically, loops are made from a contiguous list of 2d parameter space
  // curves that form a simple closed curve.
  brep.NewTrim( brep.m_E[0], false, loop, 0 ); // A to B
  brep.NewTrim( brep.m_E[1], false, loop, 1 ); // B to C
  brep.NewTrim( brep.m_E[2], false, loop, 2 ); // C to A

  // the trims are exact since we have lines on a plane
  q = brep.m_C2[0]->PointAtStart();
  //brep.m_T[0].m_P[0] = pSurface->PointAt(q.x,q.y);
  q = brep.m_C2[0]->PointAtEnd();
  //brep.m_T[0].m_P[1] = pSurface->PointAt(q.x,q.y);
  brep.m_T[0].m_type = ON_BrepTrim::boundary;
  brep.m_T[0].m_tolerance[0] = 0.0;
  brep.m_T[0].m_tolerance[1] = 0.0;

  q = brep.m_C2[0]->PointAtStart();
  //brep.m_T[1].m_P[0] = pSurface->PointAt(q.x,q.y);
  q = brep.m_C2[0]->PointAtEnd();
  //brep.m_T[1].m_P[1] = pSurface->PointAt(q.x,q.y);
  brep.m_T[1].m_type = ON_BrepTrim::boundary;
  brep.m_T[1].m_tolerance[0] = 0.0;
  brep.m_T[1].m_tolerance[1] = 0.0;

  q = brep.m_C2[0]->PointAtStart();
  //brep.m_T[2].m_P[0] = pSurface->PointAt(q.x,q.y);
  q = brep.m_C2[0]->PointAtEnd();
  //brep.m_T[2].m_P[1] = pSurface->PointAt(q.x,q.y);
  brep.m_T[2].m_type = ON_BrepTrim::boundary;
  brep.m_T[2].m_tolerance[0] = 0.0;
  brep.m_T[2].m_tolerance[1] = 0.0;

  // when debugging your code, IsValid(), IsSolid(), IsManifold() are useful
  // to check.

  model.AddModelGeometryComponent(&brep, nullptr);


  return Internal_WriteExampleModel(model, filename, error_log);
}

//int main ( int argc, const char* argv[] )
int main ()
{
    bool rc = false;
    const wchar_t* filename;
    ON::Begin();
    ON_TextLog error_log;
    ON_TextLog message_log;
    filename = L"SixEdges.3dm";
    ONX_Model* model = new ONX_Model();
    INTERNAL_INITIALIZE_MODEL(*model);

    /*
    ON_Xform trans = ON_Xform::TranslationTransformation(0, 0, -4);
    ON_Xform rotate;
    rotate.Rotation(1, 0, ON_3dVector(0, 0, 1), ON_3dPoint::Origin);
    Add_Cylinder(model, L"gold_cylinder", 2, 6, ON_Color::SaturatedGold, rotate);
    Add_Ball(model, L"blue_ball", 5, ON_Color::SaturatedBlue, trans);
    FilletSurface* fillsurface = new FilletSurface();
    ONX_ModelComponentIterator iter(*model, ON_ModelComponent::Type::ModelGeometry);
    const ON_ModelComponent* mc0 = iter.FirstComponent();
    ON_UUID uid = mc0->Id();
    ON_ModelGeometryComponent mgc = model->ModelGeometryComponentFromId(uid);
    const ON_Geometry* geo0 = mgc.Geometry(nullptr);
    const ON_NurbsSurface* nsf0 = dynamic_cast<const ON_NurbsSurface*>(geo0);
    const ON_ModelComponent* mc1 = iter.NextComponent();
    uid = mc1->Id();
    mgc = model->ModelGeometryComponentFromId(uid);
    const ON_Geometry* geo1 = mgc.Geometry(nullptr);
    const ON_NurbsSurface* nsf1 = dynamic_cast<const ON_NurbsSurface*>(geo1);
    fillsurface->SetSurface(0, *nsf0);
    fillsurface->SetSurface(1, *nsf1);
    fillsurface->SetParameter(FilletSurface::InsectCurve0Knot::v0, 0.75, FilletSurface::InsectCurve1Knot::v1, 0.7);
    fillsurface->SetRailCurve(0, 0);
    fillsurface->ReverseVector(0);
    fillsurface->ReverseVector(1);
    fillsurface->ReverseRailCurve(1);
    double knot[100] = { 0 };
    FilletSurface::SetRailUniformKnot(100, knot, knot);
    fillsurface->ComputeFillCurve(100, knot, knot);
    fillsurface->Add_FillCurves(model, L"FillCurve", ON_Color::SaturatedGreen);
    ON_NurbsSurface result_sknning;
    fillsurface->Skinning(result_sknning, 100, knot);
    const int layer_index = model->AddLayer(L"Skinning_Surface", ON_Color::SaturatedMagenta);
    model->AddManagedModelGeometryComponent(
        &result_sknning,
        Internal_CreateManagedAttributes(layer_index, L"Skinning_Surface"));
    delete fillsurface;
    */
    



    /*
    vector <ON_3dPoint> testparray;
    testparray.push_back(ON_3dPoint(0, 0, 0));
    testparray.push_back(ON_3dPoint(1, 0, 0));
    testparray.push_back(ON_3dPoint(1, 1, 0));
    testparray.push_back(ON_3dPoint(1, 1, 1));
    testparray.push_back(ON_3dPoint(0, 1.5, 1));
    testparray.push_back(ON_3dPoint(0, 2, 0));

    PiecesClothoid3d* pc3d = new PiecesClothoid3d(testparray,
        ON_3dVector(1, -0.5, 0), ON_3dVector(1, 0, -3), 8);
    pc3d->Add_to_Model(model, L"test_piece_clothoid", ON_Color::SaturatedGreen);
    */


    //Cornu_test(model);
    //EulerBsplineTest(model);
    EulerBezier2dTest(model);
    //YangMethodtest(model);

   
    rc = Internal_WriteExampleModel(*model, filename, error_log);
    if (rc)
    {
        message_log.Print(L"Successfully wrote %ls.\n", filename);
    }    
    else
    {
        message_log.Print(L"Errors while writing %ls.\n", filename);
    }  
  return 0;
}



void Add_Cylinder(ONX_Model* model, const wchar_t* name, double radius, double height, ON_Color color, ON_Xform xform)
{
    const int bIsRational = true;
    const int dim = 3;
    const int u_degree = 2;
    const int v_degree = 1;
    const int u_cv_count = 7;
    const int v_cv_count = 2;
    double u_knot[u_cv_count + u_degree - 1] = { 0,0,1.0 / 3,1.0 / 3,2.0 / 3,2.0 / 3,1,1 };
    double v_knot[v_cv_count + v_degree - 1] = { 0,1 };
    ON_3dPoint CV[u_cv_count][v_cv_count];
    for (int k = -1; k < 2; k = k + 2)
    {
        int y = (k + 1) / 2;
        CV[0][y].Set(0, -radius, k * height / 2);
        CV[1][y].Set(-radius * sqrt(3), -radius, k * height / 2);
        CV[2][y].Set(-radius * sqrt(3) / 2, radius / 2, k * height / 2);
        CV[3][y].Set(0, 2 * radius, k * height / 2);
        CV[4][y].Set(radius * sqrt(3) / 2, radius / 2, k * height / 2);
        CV[5][y].Set(radius * sqrt(3), -radius, k * height / 2);
        CV[6][y].Set(0, -radius, k * height / 2);
    }
    double w[7] = { 1,0.5,1,0.5,1,0.5,1 };

    ON_NurbsSurface nurbs_surface(dim, bIsRational,
        u_degree + 1, v_degree + 1,
        u_cv_count, v_cv_count);
    int i, j;
    for (i = 0; i < nurbs_surface.KnotCount(0); i++)
    {
        nurbs_surface.SetKnot(0, i, u_knot[i]);
    }

    for (j = 0; j < nurbs_surface.KnotCount(1); j++)
    {
        nurbs_surface.SetKnot(1, j, v_knot[j]);
    }

    for (i = 0; i < nurbs_surface.CVCount(0); i++)
    {
        for (j = 0; j < nurbs_surface.CVCount(1); j++)
        {
            nurbs_surface.SetCV(i, j, CV[i][j] * w[i]);
            nurbs_surface.SetWeight(i, j, w[i]);
        }
    }
    nurbs_surface.Transform(xform);
    ON_NurbsSurface** nubs = new ON_NurbsSurface * ();
    *nubs = new ON_NurbsSurface(nurbs_surface);
    model->AddDefaultLayer(nullptr, ON_Color::UnsetColor);
    const int layer_index = model->AddLayer(name, color);
    model->AddManagedModelGeometryComponent(
        *nubs,
        Internal_CreateManagedAttributes(layer_index, name));
    delete nubs;
}

void Add_Circle(ONX_Model* model, const wchar_t* name, double radius, ON_Color color, ON_Xform xform)
{
    ON_NurbsCurve* wiggle = new ON_NurbsCurve(
        3, // dimension
        true, // true if rational
        3,     // order = degree+1
        7      // number of control vertices
    );
    int i;

    ON_3dPoint CV[7];
    CV[0].Set(0, -radius, 0);
    CV[1].Set(-radius * sqrt(3), -radius, 0);
    CV[2].Set(-radius * sqrt(3) / 2, radius / 2, 0);
    CV[3].Set(0, 2 * radius, 0);
    CV[4].Set(radius * sqrt(3) / 2, radius / 2, 0);
    CV[5].Set(radius * sqrt(3), -radius, 0);
    CV[6].Set(0, -radius, 0);
    double weight[7] = { 1,0.5,1,0.5,1,0.5,1 };

    for (i = 0; i < 7; i++)
    {
        wiggle->SetCV(i, CV[i] * weight[i]);
        wiggle->SetWeight(i, weight[i]);
    }

    wiggle->SetKnot(0, 0.0);
    wiggle->SetKnot(1, 0.0);
    wiggle->SetKnot(2, 1.0 / 3);
    wiggle->SetKnot(3, 1.0 / 3);
    wiggle->SetKnot(4, 2.0 / 3);
    wiggle->SetKnot(5, 2.0 / 3);
    wiggle->SetKnot(6, 1.0);
    wiggle->SetKnot(7, 1.0);

    //model->AddManagedModelGeometryComponent(wiggle, nullptr);
    wiggle->Transform(xform);
    ON_NurbsCurve** nubs = new ON_NurbsCurve * ();
    *nubs = new ON_NurbsCurve(*wiggle);
    model->AddDefaultLayer(nullptr, ON_Color::UnsetColor);
    const int layer_index = model->AddLayer(name, color);
    model->AddManagedModelGeometryComponent(
        *nubs,
        Internal_CreateManagedAttributes(layer_index, name));
    delete nubs;
    delete wiggle;
}

void Add_Ball(ONX_Model* model, const wchar_t* name, double radius, ON_Color color, ON_Xform xform)
{
    const int bIsRational = true;
    const int dim = 3;
    const int u_degree = 2;
    const int v_degree = 2;
    const int u_cv_count = 7;
    const int v_cv_count = 4;
    double u_knot[u_cv_count + u_degree - 1] = { 0,0,0.25,0.5,0.5,0.75,1,1 };
    double v_knot[v_cv_count + v_degree - 1] = { 0, 0, 0.5, 1, 1 };
    ON_3dPoint CV[u_cv_count][v_cv_count];
    double weight[u_cv_count][v_cv_count];
    weight[0][0] = 1; weight[0][1] = 0.5; weight[0][2] = 0.5; weight[0][3] = 1;
  
    int i, j;
    for (i = 0; i < u_cv_count; i++)
    {
        if (i == 0 || i == 3 || i == 6)
        {
            for (j = 0; j < v_cv_count; j++)
            {
                weight[i][j] = weight[0][j];
            }
        }
        else
        {
            for (j = 0; j < v_cv_count; j++)
            {
                weight[i][j] = 0.5 * weight[0][j];
            }
        }
    }
    for (i = 0; i < u_cv_count; i++)
    {
        CV[i][0] = ON_3dPoint(0, 0, -radius) * weight[i][0];
        CV[i][3] = ON_3dPoint(0, 0, radius) * weight[i][3];
    }
    CV[0][1] = ON_3dPoint(radius, 0, -radius) * weight[0][1];
    CV[1][1] = ON_3dPoint(radius, radius, -radius) * weight[1][1];
    CV[2][1] = ON_3dPoint(-radius, radius, -radius) * weight[2][1];
    CV[3][1] = ON_3dPoint(-radius, 0, -radius) * weight[3][1];
    CV[4][1] = ON_3dPoint(-radius, -radius, -radius) * weight[4][1];
    CV[5][1] = ON_3dPoint(radius, -radius, -radius) * weight[5][1];
    CV[6][1] = ON_3dPoint(radius, 0, -radius) * weight[6][1];

    CV[0][2] = ON_3dPoint(radius, 0, radius) * weight[0][2];
    CV[1][2] = ON_3dPoint(radius, radius, radius) * weight[1][2];
    CV[2][2] = ON_3dPoint(-radius, radius, radius) * weight[2][2];
    CV[3][2] = ON_3dPoint(-radius, 0, radius) * weight[3][2];
    CV[4][2] = ON_3dPoint(-radius, -radius, radius) * weight[4][2];
    CV[5][2] = ON_3dPoint(radius, -radius, radius) * weight[5][2];
    CV[6][2] = ON_3dPoint(radius, 0, radius) * weight[6][2];

    ON_NurbsSurface nurbs_surface(dim, bIsRational,
        u_degree + 1, v_degree + 1,
        u_cv_count, v_cv_count);
    for (i = 0; i < nurbs_surface.KnotCount(0); i++)
    {
        nurbs_surface.SetKnot(0, i, u_knot[i]);
    }

    for (j = 0; j < nurbs_surface.KnotCount(1); j++)
    {
        nurbs_surface.SetKnot(1, j, v_knot[j]);
    }

    for (i = 0; i < nurbs_surface.CVCount(0); i++)
    {
        for (j = 0; j < nurbs_surface.CVCount(1); j++)
        {
            nurbs_surface.SetCV(i, j, CV[i][j]);
            nurbs_surface.SetWeight(i, j, weight[i][j]);
        }
    }
    nurbs_surface.Transform(xform);
    ON_NurbsSurface** nubs = new ON_NurbsSurface * ();
    *nubs = new ON_NurbsSurface(nurbs_surface);
    model->AddDefaultLayer(nullptr, ON_Color::UnsetColor);
    const int layer_index = model->AddLayer(name, color);
    model->AddManagedModelGeometryComponent(
        *nubs,
        Internal_CreateManagedAttributes(layer_index, name));
    delete nubs;
}

void Add_BezierCurve(ONX_Model* model, const wchar_t* name, const ON_BezierCurve* bc, ON_Color color, ON_Xform xform)
{
    ON_NurbsCurve** onc = new ON_NurbsCurve * ();
    (*onc) = new ON_NurbsCurve();
    bc->GetNurbForm(**onc);
    (*onc)->Transform(xform);
    const int layer_index = model->AddLayer(name, color);
    model->AddManagedModelGeometryComponent(
        *onc,
        Internal_CreateManagedAttributes(layer_index, name));
    delete onc;
}

void Get_Curvature_and_Torsion_of_NurbsCurve(const ON_NurbsCurve& onc, string filename)
{
    const double* knot = onc.Knot();
    int size = onc.KnotCount();
    double k0 = *knot;
    double kn = *(knot + size - 1);
    double t = 0;
    double kappa = 0;
    //double torsion = 0;
    ofstream ofs(filename);
    ON_3dPoint r;
    ON_3dVector r1;
    ON_3dVector r2;
    //ON_3dVector r3;
    for (int i = 0; i <= 1000; i++)
    {
        t = (kn - k0) / 1000 * i + k0;
        onc.Ev2Der(t, r, r1, r2);
        kappa = (ON_3dVector::CrossProduct(r1, r2)).Length() / (pow(r1.Length(), 3));
        ofs << t << " " << kappa << endl;
    }
    ofs.close();
}

void Cornu_test(ONX_Model* model)
{
    
    Cornu_Spiral* cs = new Cornu_Spiral(ON_2dPoint(0, 0), ON_2dPoint(10, 0),
        ON_2dVector(cos(10.0 / 180.0 * PI), -sin(10.0 / 180.0 * PI)),
        ON_2dVector(cos(15.0 / 180.0 * PI), sin(15.0 / 180.0 * PI))
    );
    //cs->Add_to_Model(model, L"test_Cornu_spiral", ON_Color::SaturatedGreen);
    cs->Add_Nurbs_to_Model(model, L"C_Shape", ON_Color::SaturatedMagenta);
    cs->Raise_to_3D(2, model, L"C_Shape_Rise2");
    cs->Raise_to_3D(-6, model, L"C_Shape_Rise-6");
    cs->Raise_to_3D(10, model, L"C_Shape_Rise10");
    delete cs;
    cs = new Cornu_Spiral(ON_2dPoint(0, 0), ON_2dPoint(10, 0),
        ON_2dVector(-cos(75.0 / 180.0 * PI), sin(75.0 / 180.0 * PI)),
        ON_2dVector(-cos(15.0 / 180.0 * PI), sin(15.0 / 180.0 * PI))
    );
    cs->Add_Nurbs_to_Model(model, L"S_Shape", ON_Color::SaturatedGold);
    cs->Raise_to_3D(2, model, L"S_Shape_Rise2");
    cs->Raise_to_3D(-6, model, L"S_Shape_Rise-6");
    cs->Raise_to_3D(10, model, L"S_Shape_Rise10");
    delete cs;
}

bool EulerBsplineSpiralCheck(const ON_NurbsCurve& onc)
{
    int n = onc.CVCount();
    if (n < 4)
    {
        return false;
    }
    //Compute Length
    std::vector <double> Length;
    ON_3dPoint p1, p2;
    for (int i = 0; i < n - 1; i++)
    {
        onc.GetCV(i, p1);
        onc.GetCV(i + 1, p2);
        Length.push_back((p1 - p2).Length());
    }
    //Check length is equal or not
    std::vector <double> sort_length = Length;
    std::sort(sort_length.begin(), sort_length.end());
    for (int i = 0; i < n - 2; i++)
    {
        if (abs(sort_length[i] - sort_length[i + 1]) > 1e-3)
        {
            return false;
        }
    }
    //Compute Angles
    std::vector <double> Angle;
    Angle.push_back(0);
    for (int i = 0; i < n - 2; i++)
    {
        onc.GetCV(i, p1);
        onc.GetCV(i + 1, p2);
        ON_3dVector v1 = p2 - p1;
        onc.GetCV(i + 2, p1);
        ON_3dVector v2 = p1 - p2;
        double theta = ON_3dVector::Angle(v1, v2);
        if (theta > PI / 2)
        {
            return false;
        }
        if (v1.x * v2.y - v1.y * v2.x < 0)
        {
            theta = -theta;
        }
        Angle.push_back(theta);
    }
    Angle.push_back(0);
    //Compute theta_DD
    for (int i = 2; i < n - 2; i++)
    {
        if (abs(Angle[i] * 2 - Angle[i - 1] - Angle[i + 1]) > 1e-6)
        {
            return false;
        }
    }
    //Compute delta_theta,s0,s1
    int m = n - 1;
    double Deltatheta = (Angle[m - 1] - Angle[1]) / (m - 2);
    double s0 = (m + 1) * sin(Angle[1]) + (m - 2) * sin(Angle[1] + Angle[2]) - 3 * (m - 1) * sin(Angle[1]) * cos(Angle[1]);
    double s1 = -(m + 1) * sin(Angle[n - 2]) - (m - 2) * sin(Angle[n - 3] + Angle[n - 2]) + 3 * (m - 1) * sin(Angle[n - 2]) * cos(Angle[n - 2]);
    if (Deltatheta * s0 >= 0 && s0 * s1 >= 0)
    {
        return true;
    }
    return false;
}

void SmoothingBsplineControlPolygon(ON_NurbsCurve& onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb)
{
    int n = onc.CVCount();
    int m = onc.CVCount() - 1;
    if (m <= 3)
    {
        return;
    }

    //ComputeAngle();
    std::vector <double> Angle;
    ON_3dPoint p1, p2;
    Angle.push_back(0);
    for (int i = 0; i < n - 2; i++)
    {
        onc.GetCV(i, p1);
        onc.GetCV(i + 1, p2);
        ON_3dVector v1 = p2 - p1;
        onc.GetCV(i + 2, p1);
        ON_3dVector v2 = p1 - p2;
        double theta = ON_3dVector::Angle(v1, v2);
        if (v1.x * v2.y - v1.y * v2.x < 0)
        {
            theta = -theta;
        }
        Angle.push_back(theta);
    }
    Angle.push_back(0);

    int max_count = 10000;
    for (int i = 0; i < n; i++)
    {
        if (Angle[i] > PI / 2 || Angle[i] < -PI / 2)
        {
            max_count = 2;
        }
    }
    int s_count = 1;
    double thetaDD = 10.0;
    double lengthavg = 0.0;
    onc.GetCV(0, p1);
    onc.GetCV(m, p2);
    double lengthbound = 2 * (p1 - p2).Length();

    while (s_count < max_count
        && thetaDD>1e-6
        && lengthavg < lengthbound)
    {
        vector <double> new_angle = Angle;

        for (int i = 2; i < m - 1; i++)
        {
            new_angle[i] = (new_angle[i - 1] + new_angle[i] + new_angle[i + 1]) / 3;
            onc.GetCV(i - 1, p1);
            onc.GetCV(i + 1, p2);
            ON_3dPoint re = (p2 - p1) / 2;
            re.Set(-re.y, re.x, 0);
            re = (p1 + p2) / 2 -
                re * tan(new_angle[i] / 2);
            onc.SetCV(i, re);
        }

        onc.GetCV(2, p2);
        Ta.Unitize();
        Tb.Unitize();
        double la = ON_3dVector::DotProduct(p2 - Pa, Ta);
        if (la > 0)
        {
            onc.SetCV(0, p2 - 2 * la * Ta);
        }
        else
        {
            onc.SetCV(0, Pa + la * Ta);
        }
        ON_3dVector Na = Ta;
        Na.Rotate(1, 0, ON_3dVector(0, 0, 1));
        onc.GetCV(0, p1);
        double ha = 0.5 * ON_3dVector::DotProduct(Na, Pa - 0.5 * (p1 + p2));
        onc.SetCV(1, Pa + ha * Na);

        onc.GetCV(m - 2, p2);
        double lb = ON_3dVector::DotProduct(p2 - Pb, Tb);
        if (lb < 0)
        {
            onc.SetCV(m, p2 - 2 * lb * Tb);
        }
        else
        {
            onc.SetCV(m, Pb + lb * Tb);
        }
        ON_3dVector Nb = Tb;
        Nb.Rotate(1, 0, ON_3dVector(0, 0, 1));
        onc.GetCV(m, p1);
        double hb = 0.5 * ON_3dVector::DotProduct(Nb, Pb - 0.5 * (p1 + p2));
        onc.SetCV(m - 1, Pb + hb * Nb);

        //ComputeAngle() again!
        Angle.clear();
        Angle.push_back(0);
        for (int i = 0; i < n - 2; i++)
        {
            onc.GetCV(i, p1);
            onc.GetCV(i + 1, p2);
            ON_3dVector v1 = p2 - p1;
            onc.GetCV(i + 2, p1);
            ON_3dVector v2 = p1 - p2;
            double theta = ON_3dVector::Angle(v1, v2);
            if (v1.x * v2.y - v1.y * v2.x < 0)
            {
                theta = -theta;
            }
            Angle.push_back(theta);
        }
        Angle.push_back(0);

        thetaDD = 0;
        if (m > 3)
        {
            for (int i = 2; i < m - 1; i++)
            {
                thetaDD = max(thetaDD, abs(Angle[i] * 2 - Angle[i - 1] - Angle[i + 1]));
            }
        }
        s_count += 1;
    }
}

void EulerBsplineSpiralInterpolation(ON_NurbsCurve& onc, int max_vtx_num)
{
    int v_num = onc.CVCount();
    while (v_num <= max_vtx_num && !EulerBsplineSpiralCheck(onc))
    {
        v_num = onc.CVCount();
        ON_3dPoint p1, p2;
        vector <ON_3dPoint> vp;
        onc.GetCV(0, p1);
        vp.push_back(p1);
        for (int i = 1; i < v_num; i++)
        {
            onc.GetCV(i - 1, p1);
            onc.GetCV(i, p2);
            vp.push_back(p1 * (double(i) / double(v_num)) + p2 * (1 - double(i) / double(v_num)));
        }
        onc.GetCV(v_num - 1, p1);
        vp.push_back(p1);
        double u0, u1;
        onc.GetDomain(&u0, &u1);
        onc.EvPoint(u0, p1);
        onc.EvPoint(u1, p2);
        ON_3dVector v1 =  onc.TangentAt(u0);
        ON_3dVector v2 = onc.TangentAt(u1);
        onc.Create(2, false, onc.Order(), v_num + 1);

        int i = 0;
        for (const auto &it : vp)
        {
            onc.SetCV(i, it);
            i++;
        }
        SmoothingBsplineControlPolygon(onc, p1, p2, v1, v2);
        //onc.MakeClampedUniformKnotVector();
        /*
        const double* kn = onc.Knot();
        vector <double> seeknot;
        for (int i = 0; i < onc.KnotCount(); i++)
        {
            seeknot.push_back(*(kn + i));
        }
        */
        for (int i = 0; i < onc.KnotCount(); i++)
        {
            onc.SetKnot(i, i + 1);
        }
    }
}

void EulerBsplineTest(ONX_Model* model)
{
    ON_NurbsCurve* onc = new ON_NurbsCurve();
    onc->Create(2, false, 4, 5);
    for (int i = 0; i < 7; i++)
    {
        onc->SetKnot(i, i + 1);
    }
    onc->SetCV(0, ON_3dPoint(0, 0, 0));
    onc->SetCV(1, ON_3dPoint(0, 1, 0));
    onc->SetCV(2, ON_3dPoint(1, 2, 0));
    onc->SetCV(3, ON_3dPoint(2, 4, 0));
    onc->SetCV(4, ON_3dPoint(5, 3, 0));

    const int layer_index = model->AddLayer(L"test_input_Bspline", ON_Color::SaturatedBlue);
    model->AddManagedModelGeometryComponent(
        onc,
        Internal_CreateManagedAttributes(layer_index, L"test_input_Bspline"));
    
    ON_NurbsCurve* another_onc = new ON_NurbsCurve(*onc);
    EulerBsplineSpiralInterpolation(*another_onc);
    const int layer_index2 = model->AddLayer(L"test_spiral_Bspline", ON_Color::SaturatedMagenta);
    model->AddManagedModelGeometryComponent(
        another_onc,
        Internal_CreateManagedAttributes(layer_index2, L"test_spiral_Bspline"));
        
}

void Smoothing3DControlPolygon(ON_NurbsCurve& onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb)
{
    int n = onc.CVCount();
    int m = onc.CVCount() - 1;
    if (m <= 3)
    {
        return;
    }

    //ComputeAngle();
    std::vector <double> Angle;
    ON_3dPoint p1, p2, p3;
    Angle.push_back(0);
    for (int i = 0; i < n - 2; i++)
    {
        onc.GetCV(i, p1);
        onc.GetCV(i + 1, p2);
        ON_3dVector v1 = p2 - p1;
        onc.GetCV(i + 2, p1);
        ON_3dVector v2 = p1 - p2;
        double theta = ON_3dVector::Angle(v1, v2);
        Angle.push_back(theta);
    }
    Angle.push_back(0);

    int max_count = 10000;
    for (int i = 0; i < n; i++)
    {
        if (Angle[i] > PI / 2 || Angle[i] < -PI / 2)
        {
            max_count = 2;
        }
    }
    int s_count = 1;
    double thetaDD = 10.0;
    double lengthavg = 0.0;
    onc.GetCV(0, p1);
    onc.GetCV(m, p2);
    double lengthbound = 2 * (p1 - p2).Length();

    //Compute Normal Vector
    std::vector <ON_3dVector> normal;
    normal.push_back(ON_3dVector(0, 0, 0));
    for (int i = 0; i < n - 2; i++)
    {
        onc.GetCV(i, p1);
        onc.GetCV(i + 1, p2);
        onc.GetCV(i + 2, p3);
        ON_3dVector vv = ON_3dVector::CrossProduct(p2 - p1, p3 - p2);
        //vv.Unitize();
        normal.push_back(vv);
    }
    normal.push_back(ON_3dVector(0, 0, 0));


    while (s_count < max_count
        && thetaDD>1e-6
        && lengthavg < lengthbound)
    {
        //Compute New Vertex 2--n-2
        vector <double> new_angle = Angle;
        vector <ON_3dVector> new_normal = normal;

        for (int i = 2; i < m - 1; i++)
        {
            new_angle[i] = (new_angle[i - 1] + new_angle[i] + new_angle[i + 1]) / 3;
            onc.GetCV(i - 1, p1);
            onc.GetCV(i + 1, p2);
            new_normal[i] = (new_normal[i - 1] + new_normal[i] + new_normal[i + 1]) / 3;
            ON_3dVector beta = ON_3dVector::CrossProduct(p2 - p1, new_normal[i]);
            beta.Unitize();
            ON_3dPoint re = (p1 + p2) / 2 + beta * tan(new_angle[i] / 2);
            onc.SetCV(i, re);
        }



    }
}