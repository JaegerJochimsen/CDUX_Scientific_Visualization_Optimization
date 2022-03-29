#include <iostream>
#include "util.h"

void foo(){if(SUBSET) std::cout << "imported\n";}

int vectorIndex(std::vector<Vertex> v, Vertex p){
    std::vector<Vertex>::iterator itr = std::find(v.begin(), v.end(), p);
    return std::distance(v.begin(), itr);
}

/* pointInBoundingBox()
 * 
 * Description:
 *      This function determines if a given point is within one of the bounding boxes defined by BXMin, BXMax, BYMin, BYMax, ...
 *
 * Return:
 *      true if the point in the box, false otherwise
 */

bool pointInBoundingBox(double vx, double vy, double vz, float bxMin, float bxMax, float byMin, float byMax, float bzMin, float bzMax){
	bool contains = (double)bxMin <= vx && vx <= (double)bxMax &&
					(double)byMin <= vy && vy <= (double)byMax &&
					(double)bzMin <= vz && vz <= (double)bzMax;
	return contains;
}

int getDomainIndex(double x, double y, double z){
	if(!pointInBoundingBox(x,y,z,BXMin, BXMax, BYMin, BYMax, BZMin, BZMax)) return -1;
    return (int)(x/W) + (int)(y/W)*16 + (int)(z/W)*4;
}

std::vector<int> neighborDoms(Vertex v){
    // Place point in all adjacent subdomains and then find the subDomain ID for that
    // This will remain -1 if it is out of bounds
    std::vector<int> neighborDomains(26, -1);

    neighborDomains[0] = getDomainIndex(v.x(), v.y(), v.z() +W);
    neighborDomains[1] = getDomainIndex(v.x(), v.y(), v.z() -W);
    neighborDomains[2] = getDomainIndex(v.x(), v.y() +W, v.z());
    neighborDomains[3] = getDomainIndex(v.x(), v.y() +W, v.z() +W);
    neighborDomains[4] = getDomainIndex(v.x(), v.y() +W, v.z() -W);
    neighborDomains[5] = getDomainIndex(v.x(), v.y() -W, v.z());
    neighborDomains[6] = getDomainIndex(v.x(), v.y() -W, v.z() +W);
    neighborDomains[7] = getDomainIndex(v.x(), v.y() -W, v.z() -W);

    neighborDomains[8] = getDomainIndex(v.x() +W, v.y(), v.z());
    neighborDomains[9] = getDomainIndex(v.x() +W, v.y(), v.z() +W);
    neighborDomains[10] = getDomainIndex(v.x() +W, v.y(), v.z() -W);
    neighborDomains[11] = getDomainIndex(v.x() +W, v.y() +W, v.z());
    neighborDomains[12] = getDomainIndex(v.x() +W, v.y() +W, v.z() +W);
    neighborDomains[13] = getDomainIndex(v.x() +W, v.y() +W, v.z() -W);
    neighborDomains[14] = getDomainIndex(v.x() +W, v.y() -W, v.z());
    neighborDomains[15] = getDomainIndex(v.x() +W, v.y() -W, v.z() +W);
    neighborDomains[16] = getDomainIndex(v.x() +W, v.y() -W, v.z() -W);

    neighborDomains[17] = getDomainIndex(v.x() -W, v.y(), v.z());
    neighborDomains[18] = getDomainIndex(v.x() -W, v.y(), v.z() +W);
    neighborDomains[19] = getDomainIndex(v.x() -W, v.y(), v.z() -W);
    neighborDomains[20] = getDomainIndex(v.x() -W, v.y() +W, v.z());
    neighborDomains[21] = getDomainIndex(v.x() -W, v.y() +W, v.z() +W);
    neighborDomains[22] = getDomainIndex(v.x() -W, v.y() +W, v.z() -W);
    neighborDomains[23] = getDomainIndex(v.x() -W, v.y() -W, v.z());
    neighborDomains[24] = getDomainIndex(v.x() -W, v.y() -W, v.z() +W);
    neighborDomains[25] = getDomainIndex(v.x() -W, v.y() -W, v.z() -W);

    return neighborDomains;
}

double euclidDistance3D(double *p0, double *p1){
    return sqrt(pow(p0[0] - p1[0], 2) + pow(p0[1] - p1[1], 2) + pow(p0[2] - p1[2], 2));
}

/*
 * This is based off a paper titled:  "3D Barycentric Interpolation Method for Animal Brainmapping"
 *
 * Citation:
 *       P. Pešková, M. Piorecký, Č. Vejmola, V. Koudelka, V. Krajča and V. Piorecká, "3D Barycentric Interpolation Method for Animal Brainmapping," 2019
 *       E-Health and Bioengineering Conference (EHB), 2019, pp. 1-4, doi: 10.1109/EHB47216.2019.8969996.
 * 
 * Interpolation
 *      
 *              Σ (V_e/d_ie)
 *      V_i =  -------------  
 *              Σ (1/d_ie)
 *
 *      V_i is the interpolated vertex
 *      V_e is (in turn) each of the 4 vertices of the tetrahedron that contains V_i's start vert
 *      d_ie is the distance from the current base vertex V_e and V_i
 */

void barycentricInterp(double *endPt, double *pt, 
                       double *p0Start, double *p1Start, double *p2Start, double *p3Start, 
                       double *p0End, double *p1End, double *p2End, double *p3End)
{
    // determine distances for barycentric interp
    double d0new = euclidDistance3D(pt, p0Start);
    double d1new = euclidDistance3D(pt, p1Start);
    double d2new = euclidDistance3D(pt, p2Start);
    double d3new = euclidDistance3D(pt, p3Start);
    double total = 1.0F/d0new + 1.0F/d1new + 1.0F/d2new + 1.0F/d3new;
    
    // order 2 (1) interpolation seems to split the error better than order 3 (in my limited testing)
    double t0Prop = d0new*d0new;
    double t1Prop = d1new*d1new;
    double t2Prop = d2new*d2new;
    double t3Prop = d3new*d3new;

    // Σ (V_e/d_ie)
    double numerator0[3] = {p0End[0]/t0Prop, p0End[1]/t0Prop, p0End[2]/t0Prop};     
    double numerator1[3] = {p1End[0]/t1Prop, p1End[1]/t1Prop, p1End[2]/t1Prop};
    double numerator2[3] = {p2End[0]/t2Prop, p2End[1]/t2Prop, p2End[2]/t2Prop};
    double numerator3[3] = {p3End[0]/t3Prop, p3End[1]/t3Prop, p3End[2]/t3Prop};
    
    // Σ (1/d_ie)
    double lowerSum = 1.0/t0Prop + 1.0/t1Prop + 1.0/t2Prop + 1.0/t3Prop; 

    // [Σ (V_e/d_ie)] / [Σ (1/d_ie)]
    double newX = (numerator0[0] + numerator1[0] + numerator2[0] + numerator3[0])/lowerSum;
    double newY = (numerator0[1] + numerator1[1] + numerator2[1] + numerator3[1])/lowerSum;
    double newZ = (numerator0[2] + numerator1[2] + numerator2[2] + numerator3[2])/lowerSum;

    // store new interpolated value
    endPt[0] = newX;
    endPt[1] = newY;
    endPt[2] = newZ;

    return;
}

void sortInputs(vtkDataSetReader *               inputRdr, 
                std::vector<std::vector<Vertex>> &startBasis, 
                std::vector<std::vector<Vertex>> &endBasis,
                std::vector<std::vector<Vertex>> &queryPts,
                std::vector<Vertex>              &allStart,
                long                             numBasisPts)
{
    // a bounding box per subdomain
    int ***start_bboxes = (int ***)malloc(sizeof(int **)*64);
    for(int j = 0; j < 64; ++j){
        start_bboxes[j] = (int **)malloc(sizeof(int *)*8);
        for(int k = 0; k < 8; ++k){
            start_bboxes[j][k] = (int *)calloc(3, sizeof(int));
        } 
    }

    int subdomain_count = 0;
    long i = 0;
    for(; i < (long)(numBasisPts/2); ++i){
        double *start_point = inputRdr->GetOutput()->GetPoint(2*i);
        double *end_point = inputRdr->GetOutput()->GetPoint(2*i + 1);
 
        // get subdomain id 
        int subdomain_id = getDomainIndex(start_point[0], start_point[1], start_point[2]);

        // Take random subset of input pts (of size SUBSET)
        //double r = ((double)rand() / (RAND_MAX));
        //double r = 0.0;
        if(subdomain_count < SUBSET){
            // query pts are just the start, their projected end should be the same for this test
            queryPts[subdomain_id].push_back(Vertex(start_point[0], start_point[1], start_point[2]));
            subdomain_count++;
        }
        else{
            // for use in creating global delaunay
            allStart.push_back(Vertex(start_point[0], start_point[1], start_point[2]));
            
            // for use in creating local and neighbor delaunays
            startBasis[subdomain_id].push_back(Vertex(start_point[0], start_point[1], start_point[2]));
            endBasis[subdomain_id].push_back(Vertex(end_point[0], end_point[1], end_point[2]));
        }
    }
    return;
}
