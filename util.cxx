#include <iostream>
#include "util.h"

int vectorIndex(std::vector<Vertex> v, Vertex p){
    std::vector<Vertex>::iterator itr = std::find(v.begin(), v.end(), p);
    return std::distance(v.begin(), itr);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

bool pointInBoundingBox(double vx, double vy, double vz, double bxMin, 
                        double bxMax, double byMin, double byMax, 
                        double bzMin, double bzMax){

	bool contains = (double)bxMin <= vx && vx <= (double)bxMax &&
					(double)byMin <= vy && vy <= (double)byMax &&
					(double)bzMin <= vz && vz <= (double)bzMax;
	return contains;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int getDomainIndex(double x, double y, double z){
	if(!pointInBoundingBox(x,y,z,BXMin, BXMax, BYMin, BYMax, BZMin, BZMax)) 
        return -1;
    return (int)(x/W) + (int)(y/W)*16 + (int)(z/W)*4;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

std::vector<int> neighborDoms(Vertex v){
    // Place point in all adjacent subdomains and then find the subDomain ID 
    // for that.
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

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

double euclidDistance3D(double *p0, double *p1){
    return sqrt(pow(p0[0] - p1[0], 2) + 
                pow(p0[1] - p1[1], 2) + 
                pow(p0[2] - p1[2], 2));
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void barycentricInterp(double *endPt, double *pt, 
                       double *p0Start, double *p1Start, 
                       double *p2Start, double *p3Start, 
                       double *p0End, double *p1End, 
                       double *p2End, double *p3End)
{
    // determine distances for barycentric interp
    double d0new = euclidDistance3D(pt, p0Start);
    double d1new = euclidDistance3D(pt, p1Start);
    double d2new = euclidDistance3D(pt, p2Start);
    double d3new = euclidDistance3D(pt, p3Start);
    double total = 1.0F/d0new + 1.0F/d1new + 1.0F/d2new + 1.0F/d3new;
    
    // order 2 interpolation seems to split the error better than 
    // order 1 or 3 (in my limited testing)
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

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void modifyBoundingBox(double *bbox, double *pt){
    // if pt is a new smallest
    if(bbox[0] > pt[0])
        bbox[0] = pt[0];

    if(bbox[1] > pt[1])
        bbox[1] = pt[1];

    if(bbox[2] > pt[2])
        bbox[2] = pt[2];

    // if pt is new biggest
    if(bbox[3] < pt[0])
        bbox[3] = pt[0];

    if(bbox[4] < pt[1])
        bbox[4] = pt[1];

    if(bbox[5] < pt[2])
        bbox[5] = pt[2];
    return;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

bool edgePoint(double vx, double vy, double vz){
    
    // get the bottom left vertex
    double bottom_x = (int)(vx/W)*W;
    double bottom_y = (int)(vy/W)*W;
    double bottom_z = (int)(vz/W)*W;
    
    double outerbbox[] = { bottom_x, bottom_y, bottom_z,
                           bottom_x + (double)W, bottom_y + (double)W, bottom_z + (double)W,
                        };

    double center[] = {bottom_x + (double)(W/2), bottom_y + (double)(W/2), bottom_z + (double)(W/2)};

    double innerbbox[] = {center[0] - (double)N*(double)(W/2), center[1] - (double)N*(double)(W/2), center[2] - (double)N*(double)(W/2),
                          center[0] + (double)N*(double)(W/2), center[1] + (double)N*(double)(W/2), center[2] + (double)N*(double)(W/2), 
                        };
    
    // if in outer box and not in inner box
    bool in_outer = pointInBoundingBox(vx,vy,vz,outerbbox[0],outerbbox[3],outerbbox[1], outerbbox[4], outerbbox[2],outerbbox[5]);
    bool in_inner = pointInBoundingBox(vx,vy,vz,innerbbox[0],innerbbox[3], innerbbox[1], innerbbox[4], innerbbox[2],innerbbox[5]);

    return in_outer && !in_inner;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int *nearestSubdomains(double vx, double vy, double vz){
    int *nearest = (int *)calloc((size_t)26, sizeof(int));
    for(int i = 0; i < 26; ++i) nearest[i] = -1;
    
    // Here we will define a new convention for face naming
    // These will be used to fill values for bounding boxes (1 per face)
    // Face # - Orientation
    // 0        Bottom
    // 1        Left
    // 2        Front
    // 3        Right
    // 4        Back
    // 5        Top

    // copy globals as an appropriate type
    double w = (double)W;
    double n = (double)N;
    double comp_n = (double)(1.0 - n); // compliment of N

    // get the bottom left vertex of subdomain
    double bottom_x = (int)(vx/W)*W;
    double bottom_y = (int)(vy/W)*W;
    double bottom_z = (int)(vz/W)*W;
    
    // Define bounding boxes. 
    // Ordered coord0, coordmax0, ...                  
    double bottom[] = {bottom_x, bottom_x + w,         
                       bottom_y, bottom_y + comp_n*w,  
                       bottom_z, bottom_z + w
                      };     

    double left[] = {bottom_x, bottom_x + comp_n*w,
                     bottom_y, bottom_y + w,
                     bottom_z, bottom_z + w
                    };

    double front[] = {bottom_x, bottom_x + w,
                      bottom_y, bottom_y + w,
                      bottom_z, bottom_z + comp_n*w
                     };

    double right[] = {bottom_x + n*w, bottom_x + w,
                      bottom_y, bottom_y + w,
                      bottom_z, bottom_z + w
                     };

    double back[] = {bottom_x, bottom_x + w,
                     bottom_y, bottom_y + w,
                     bottom_z + n*w, bottom_z + w 
                    };

    double top[] = {bottom_x, bottom_x + w,
                    bottom_y + n*w, bottom_y + w,
                    bottom_z, bottom_z + w
                   };

    int in_bottom = (int)pointInBoundingBox(vx, vy, vz, bottom[0], bottom[1], 
                                           bottom[2], bottom[3],
                                           bottom[4], bottom[5]);

    int in_left = (int)pointInBoundingBox(vx, vy, vz, left[0], left[1], 
                                           left[2], left[3],
                                           left[4], left[5]);

    int in_front = (int)pointInBoundingBox(vx, vy, vz, front[0], front[1], 
                                           front[2], front[3],
                                           front[4], front[5]);

    int in_right = (int)pointInBoundingBox(vx, vy, vz, right[0], right[1], 
                                           right[2], right[3],
                                           right[4], right[5]);

    int in_back = (int)pointInBoundingBox(vx, vy, vz, back[0], back[1], 
                                           back[2], back[3],
                                           back[4], back[5]);

    int in_top = (int)pointInBoundingBox(vx, vy, vz, top[0], top[1], 
                                                     top[2], top[3],
                                                     top[4], top[5]);

    // bottom face
    if(in_bottom) nearest[0] = getDomainIndex(vx,vy - w,vz);

    // left face
    if(in_left) nearest[1] = getDomainIndex(vx - w,vy,vz);

    // front face
    if(in_front) nearest[2] = getDomainIndex(vx,vy,vz - w);

    // right face
    if(in_right) nearest[3] = getDomainIndex(vx + w,vy,vz);

    // back face
    if(in_back) nearest[4] = getDomainIndex(vx,vy,vz + w);

    // top face
    if(in_top) nearest[5] = getDomainIndex(vx,vy + w,vz);


    // bottom left diag
    if(in_bottom && in_left) nearest[6] = getDomainIndex(vx - w, vy - w, vz);

    // bottom right diag
    if(in_bottom && in_right) nearest[7] = getDomainIndex(vx + w,vy - w,vz);

    // bottom front diag
    if(in_bottom && in_front) nearest[8] = getDomainIndex(vx,vy - w,vz - w);

    // bottom back diag
    if(in_bottom && in_back) nearest[9] = getDomainIndex(vx,vy - w,vz + w);

    // top left diag
    if(in_top && in_left) nearest[10] = getDomainIndex(vx - w,vy + w,vz);
    printf("past 12th\n");

    // top right diag
    if(in_top && in_right) nearest[11] = getDomainIndex(vx + w,vy + w,vz);

    // top front diag
    if(in_top && in_front) nearest[12] = getDomainIndex(vx,vy + w,vz - w);

    // top back diag
    if(in_top && in_back) nearest[13] = getDomainIndex(vx,vy + w,vz + w);

    // left back diag
    if(in_left && in_back) nearest[14] = getDomainIndex(vx - w, vy, vz + w);

    // left front diag
    if(in_left && in_front) nearest[15] = getDomainIndex(vx - w, vy, vz - w);

    // right back diag
    if(in_right && in_back) nearest[16] = getDomainIndex(vx + w, vy, vz - w);

    // right front diag
    if(in_right && in_front) nearest[17] = getDomainIndex(vx + w, vy, vz + w);


    // left front top diag
    if(in_left && in_front && in_top) nearest[18] = getDomainIndex(vx - w, vy + w, vz - w);

    // left back top diag
    if(in_left && in_back && in_top) nearest[19] = getDomainIndex(vx - w, vy + w, vz + w);

    // left back bottom diag
    if(in_left && in_back && in_bottom) nearest[20] = getDomainIndex(vx - w, vy - w, vz + w);

    // left front bottom diag
    if(in_left && in_front && in_bottom) nearest[21] = getDomainIndex(vx - w, vy - w, vz - w);
    
    // right front top diag
    if(in_right && in_front && in_top) nearest[22] = getDomainIndex(vx + w, vy + w, vz - w);

    // right back top diag
    if(in_right && in_back && in_top) nearest[23] = getDomainIndex(vx + w, vy + w, vz + w);

    // right back bottom diag
    if(in_right && in_back && in_bottom) nearest[24] = getDomainIndex(vx + w, vy - w, vz + w);

    // right front bottom diag
    if(in_right && in_front && in_bottom) nearest[25] = getDomainIndex(vx + w, vy - w, vz - w);

    return nearest;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void sortInputs(vtkDataSetReader *               inputRdr, 
                std::vector<std::vector<Vertex>> &startBasis, 
                std::vector<std::vector<Vertex>> &endBasis,
                std::vector<std::vector<Vertex>> &queryPts,
                std::vector<Vertex>              &allStart,
                long                             numBasisPts,
                bool                             doGlobal, 
                bool                             doLocal, 
                bool                             doNeighbor 
                )
{
    int subdomain_count = 0;
    long i = 0;
    for(; i < (long)(numBasisPts/2); ++i){
        double *start_point = inputRdr->GetOutput()->GetPoint(2*i);
        double *end_point = inputRdr->GetOutput()->GetPoint(2*i + 1);
 
        // get subdomain id 
        int subdomain_id = getDomainIndex(start_point[0], start_point[1], start_point[2]);
        printf("subdomain: %d\n", subdomain_id);
        if(subdomain_id < 0) continue;
        if(subdomain_id > 63) printf("HERE!\n");
        // Take subset of input pts (of size SUBSET)
        if(subdomain_count < SUBSET){
            // query pts are just the start, their projected end should be the same for this test
            queryPts[subdomain_id].push_back(Vertex(start_point[0], start_point[1], start_point[2]));
            subdomain_count++;
        }
        else{
            if(doLocal || doNeighbor){
                startBasis[subdomain_id].push_back(Vertex(start_point[0], start_point[1], start_point[2]));
                //endBasis[subdomain_id].push_back(Vertex(start_point[0], start_point[1], start_point[2]));
            }

            if(doNeighbor && edgePoint(start_point[0], start_point[1], start_point[2])){ 
                printf("after edgePoint()\n");
                int *nearest = nearestSubdomains(start_point[0], start_point[1], start_point[2]);
                printf("after nearest\n");
                for(int k = 0; k < 26; ++k){

                    if(nearest[k] == -1) continue; // since we are adding to nearest sequentially, once we see a -1 we're done

                    printf("nearest[%d]=%d\n", k, nearest[k]);

                    startBasis[nearest[k]].push_back(Vertex(start_point[0], start_point[1], start_point[2]));
                    //endBasis[nearest[k]].push_back(Vertex(start_point[0], start_point[1], start_point[2]));
                }
                free(nearest);
            }
  
            if(doGlobal){
                // for use in creating global delaunay
                allStart.push_back(Vertex(start_point[0], start_point[1], start_point[2]));
            }
        }
    }
    return;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void buildNeighborhood(std::vector<Vertex> &neighborhood, std::vector<std::vector<Vertex>> start_basis, double ***bboxes, int domId){
    // add the central subdomain
    neighborhood.insert(neighborhood.begin(), start_basis[domId].begin(), start_basis[domId].end());
    
    // check relationship between query box and central domain box
    // do so only in the direction we need to
    // Here identify which directions we need to expand in
    int lx = int(bboxes[1][domId][0] < bboxes[0][domId][0]);
    int ly = int(bboxes[1][domId][1] < bboxes[0][domId][1]);
    int lz = int(bboxes[1][domId][2] < bboxes[0][domId][2]);
    int lcase = lx + ly*2 + lz*4;
    int  id;
    std::vector<int> ids;

    switch(lcase){
        case 0:
            // point is inside box no need to expand it
            break;
        case 1:
            // need to expand in the negative x direction
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            break;
        case 2:
            // need to expand in the negative y direction
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] - W, bboxes[0][domId][2]);
            ids.push_back(id);
            break;
        case 3:
            // need to expand in the negative x AND negative y directions (and the XY diagonal)
            // X
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            // Y
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] - W, bboxes[0][domId][2]);
            ids.push_back(id);
            // XY diagonal
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1] - W, bboxes[0][domId][2]);
            ids.push_back(id);
            break;
        case 4:
            // need to expand in the negative z direction
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] - W);
            ids.push_back(id);
            break;
        case 5:
            // need to expand in the negative z AND negative x directions (and the XZ diagonal)
            // X
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            // Z
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] - W);
            ids.push_back(id);
            // XZ diagonal
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1], bboxes[0][domId][2] - W);
            ids.push_back(id);
            break;
        case 6:
            // need to expand in the negative y AND negative z directions (and along the YZ diagonal)
            // Y
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] - W, bboxes[0][domId][2]);
            ids.push_back(id);
            // Z
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] - W);
            ids.push_back(id);
            // YZ diag
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] - W, bboxes[0][domId][2] - W);
            ids.push_back(id);
            break;
        case 7:
            // expand in all dirs x, y, and z
            // X
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            // Y
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] - W, bboxes[0][domId][2]);
            ids.push_back(id);
            // Z
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] - W);
            ids.push_back(id);
            // YZ diag
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] - W, bboxes[0][domId][2] - W);
            ids.push_back(id);
            // XZ diag
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1], bboxes[0][domId][2] - W);
            ids.push_back(id);
            // XY diag
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1] - W, bboxes[0][domId][2]);
            ids.push_back(id);
            // XYZ diag
            id = getDomainIndex(bboxes[0][domId][0] - W, bboxes[0][domId][1] - W, bboxes[0][domId][2] - W);
            ids.push_back(id);
            break;
    }
    
    int ux = int(bboxes[1][domId][3] > bboxes[0][domId][3]);
    int uy = int(bboxes[1][domId][4] > bboxes[0][domId][4]);
    int uz = int(bboxes[1][domId][5] > bboxes[0][domId][5]);
    int ucase = ux + uy*2 + uz*4;
    id = -1;

    switch(ucase){
        case 0:
            // point is inside box no need to expand it
            break;
        case 1:
            // need to expand in the positive x direction
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            break;
        case 2:
            // need to expand in the positive y direction
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] + W, bboxes[0][domId][2]);
            ids.push_back(id);
            break;
        case 3:
            // need to expand in the positive x AND positive y directions (and the XY diagonal)
            // X
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            // Y
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] + W, bboxes[0][domId][2]);
            ids.push_back(id);
            // XY diagonal
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1] + W, bboxes[0][domId][2]);
            ids.push_back(id);
            break;
        case 4:
            // need to expand in the positive z direction
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] + W);
            ids.push_back(id);
            break;
        case 5:
            // need to expand in the positive z AND positive x directions (and the XZ diagonal)
            // X
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            // Z
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] + W);
            ids.push_back(id);
            // XZ diagonal
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1], bboxes[0][domId][2] + W);
            ids.push_back(id);
            break;
        case 6:
            // need to expand in the positive y AND positive z directions (and along the YZ diagonal)
            // Y
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] + W, bboxes[0][domId][2]);
            ids.push_back(id);
            // Z
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] + W);
            ids.push_back(id);
            // YZ diag
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] + W, bboxes[0][domId][2] + W);
            ids.push_back(id);
            break;
        case 7:
            // expand in all dirs x, y, and z
            // X
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1], bboxes[0][domId][2]);
            ids.push_back(id);
            // Y
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] + W, bboxes[0][domId][2]);
            ids.push_back(id);
            // Z
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1], bboxes[0][domId][2] + W);
            ids.push_back(id);
            // YZ diag
            id = getDomainIndex(bboxes[0][domId][0], bboxes[0][domId][1] + W, bboxes[0][domId][2] + W);
            ids.push_back(id);
            // XZ diag
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1], bboxes[0][domId][2] + W);
            ids.push_back(id);
            // XY diag
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1] + W, bboxes[0][domId][2]);
            ids.push_back(id);
            // XYZ diag
            id = getDomainIndex(bboxes[0][domId][0] + W, bboxes[0][domId][1] + W, bboxes[0][domId][2] + W);
            ids.push_back(id);
            break;
    }
  
    // grab all IDs for neighbor domains
    //std::vector<int> neighborhoodIDs = neighborDoms(subdomain_start_basis[domId][0]);
    //
    // add the neighboring subdomains
    //int count = 0;
    for(int i = 0; i < ids.size(); ++i){
    
        int curr_id = ids[i];
        
        if(curr_id != -1 && start_basis[curr_id].size() != 0){
                neighborhood.insert(neighborhood.begin(), start_basis[curr_id].begin(), start_basis[curr_id].end());
        }
    }
    return;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

