#ifndef _UTIL_H_
#define _UTIL_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkTetra.h>

#define SUBSET 1000
#define W 2.5
#define BXMin 0
#define BYMin 0
#define BZMin 0
#define BXMax 10
#define BYMax 10
#define BZMax 10

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K> Vb; 
typedef CGAL::Triangulation_data_structure_3<Vb, CGAL::Triangulation_cell_base_3<K>, CGAL::Parallel_tag> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Point Vertex;


void foo();

int vectorIndex(std::vector<Vertex>, Vertex);

bool pointInBoundingBox(double vx, double vy, double vz, float bxMin, float bxMax, float byMin, float byMax, float bzMin, float bzMax);

int getDomainIndex(double x, double y, double z);

std::vector<int> neighborDoms(Vertex v);

double euclidDistance3D(double *p0, double *p1);

void barycentricInterp(double *endPt, double *pt, 
                       double *p0Start, double *p1Start, double *p2Start, double *p3Start, 
                       double *p0End, double *p1End, double *p2End, double *p3End);

void sortInputs(vtkDataSetReader *               inputRdr, 
                std::vector<std::vector<Vertex>> &startBasis, 
                std::vector<std::vector<Vertex>> &endBasis,
                std::vector<std::vector<Vertex>> &queryPts,
                std::vector<Vertex>              &allStart,
                long                             numBasisPts);

#endif /* _UTIL_H_ */
