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
#define BXMin 0.0
#define BYMin 0.0
#define BZMin 0.0
#define BXMax 10.0
#define BYMax 10.0
#define BZMax 10.0
#define N 0.9

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K> Vb; 
typedef CGAL::Triangulation_data_structure_3<Vb, CGAL::Triangulation_cell_base_3<K>, CGAL::Parallel_tag> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Point Vertex;


/* vectorIndex()
 * Description: find the index of Vertex p in vector v
 * Return: index
 */
int vectorIndex(std::vector<Vertex> v, Vertex p);



/* pointInBoundingBox()
 * 
 * Description:
 *      This function determines if a given point is within one of the 
 *      bounding boxes defined by BXMin, BXMax, BYMin, BYMax, ...
 *
 * Return:
 *      true if the point in the box, false otherwise
 */
bool pointInBoundingBox(double vx, double vy, double vz, double bxMin, 
                        double bxMax, double byMin, double byMax, 
                        double bzMin, double bzMax);



/* getDomainIndex()
 * Description: get the domain index (currently 0-63) that the 
 *              point (x,y,z) is in
 * Return: index
 */
int getDomainIndex(double x, double y, double z);



/* neighborDoms()
 * Description: collect the indexes of the neighboring subdomains for a given
 *              subdomain.
 * Return: vector containing all the indexes of neighbor subdomains
 */
std::vector<int> neighborDoms(Vertex v);



/* euclidDistance3D()
 * Description: 3D distance between two points
 *
 */
double euclidDistance3D(double *p0, double *p1);



/*
 * This is based off a paper titled:  
 * "3D Barycentric Interpolation Method for Animal Brainmapping"
 *
 * Citation:
 *       P. Pešková, M. Piorecký, Č. Vejmola, V. Koudelka, V. Krajča 
 *       and V. Piorecká, 
 *       "3D Barycentric Interpolation Method for Animal Brainmapping," 2019
 *       E-Health and Bioengineering Conference (EHB), 2019, pp. 1-4, 
 *       doi: 10.1109/EHB47216.2019.8969996.
 * 
 * Interpolation
 *      
 *              Σ (V_e/d_ie)
 *      V_i =  -------------  
 *              Σ (1/d_ie)
 *
 *      V_i is the interpolated vertex
 *
 *      V_e is (in turn) each of the 4 vertices of the tetrahedron that 
 *      contains V_i's start vert
 *
 *      d_ie is the distance from the current base vertex V_e and V_i
 */
void barycentricInterp(double *endPt, double *pt, 
                       double *p0Start, double *p1Start, 
                       double *p2Start, double *p3Start, 
                       double *p0End, double *p1End, double *p2End, 
                       double *p3End);



/* modifyBoudningBox()
 * Description: expand a bounding box if a point is outside it
 * Return: void, modifies parameter bbox
 */
void modifyBoundingBox(double **bbox, double *pt);



/* sortInputs()
 * Description: sort input points into start, end, and query points and 
 *              place them in the correct vector (which corresponds to 
 *              the subdomain they are located in)
 *
 * Return: none
 * Side Effect: fills passed in vectors
 */
void sortInputs(vtkDataSetReader *               inputRdr, 
                std::vector<std::vector<Vertex>> &startBasis, 
                std::vector<std::vector<Vertex>> &endBasis,
                std::vector<std::vector<Vertex>> &queryPts,
                std::vector<Vertex>              &allStart,
                long                             numBasisPts,
                bool                             doGlobal,
                bool                             doLocal,
                bool                             doNeighbor);

void buildNeighborhood(std::vector<Vertex> &neighborhood, std::vector<std::vector<Vertex>> start_basis, double ***bboxes, int domId);


/* edgePoint()
 * 
 * Description: determine if a point is on the border of its subdomain. 
 * The distance from the border that is the "cutoff" is determined by the macro 
 * N.
 *
 * Return: true if on edge, false otherwise
 */
bool edgePoint(double vx, double vy, double vz);

/* nearestSubdomains()
 *
 * Description: determine which edges the point is closest to. Intended to be
 * used on an edge point.
 *
 * Return: an int[] with ids for the nearest subdomains
 *
 */
int *nearestSubdomains(double vx, double vy, double vz);

#endif /* _UTIL_H_ */
