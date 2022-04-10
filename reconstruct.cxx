#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <iterator> // for vector find
#include <unistd.h> // for strtok
#include <cmath>
#include <cstdlib>
#include <string>
#include <pthread.h>
#include "util.h"

#define THREADPOOL_SIZE 64
typedef struct thread_args 
{ 
    std::vector<Vertex> starts;
    std::vector<Vertex> query;
    std::vector<double> triang_times;
    std::vector<double> locate_times;
    std::vector<double> interp_times;
    std::vector<double> errors;
} Thread_Args;


int main(int argc, char *argv[]){
    // ---------------------------------------------------------------//
    // Parse flags for delaunay style (global, local, neighborhood) 
    // ---------------------------------------------------------------//
    bool parallel = false, doGlobal = false, doLocal = false, doNeighbor = false;
    int opt, optind;

    opterr = 0;
    while((opt = getopt(argc, argv, "pnlg")) != -1){
        switch(opt){
            case 'g': doGlobal = true; break;
            case 'l': doLocal = true; break;
            case 'p': parallel = true; break;
            case 'n': doNeighbor = true; break;
            default:
                      std::cout << argv[0] << ": illegal option, '-" << optopt << "'" << '\n';
                      std::cerr << "usage: " << argv[0] << " [-p] [-l] [-n] [-g] " << '\n';
                      return EXIT_FAILURE;
        }
    }

    // check flag conditions
    if(doLocal && doGlobal){
        std::cerr << "FlagError: cannot specify -g (global) and -l (local) simultaneously" << '\n';
        return EXIT_FAILURE;
    }

    if(doLocal && doNeighbor){
        std::cerr << "FlagError: cannot specify -n (neighbor) and -l (local) simultaneously" << '\n';
        return EXIT_FAILURE;
    }

    if(doNeighbor && doGlobal){
        std::cerr << "FlagError: cannot specify -g (global) and -n (neighbor) simultaneously" << '\n';
        return EXIT_FAILURE;
    }

    // default to local Delaunay
    if(!doLocal && !doGlobal && !doNeighbor){
        doLocal = true;
    }

    // ---------------------------------------------------------------//
    // Organize points into start and end (these are basis flows)
    //
    // Note: query points will be equivalent to the start points for 
    //       accuracy testing
    // ---------------------------------------------------------------//
    vtkDataSetReader *inputRdr = vtkDataSetReader::New();
    inputRdr->SetFileName("Groundtruth_0_0.vtk");
    inputRdr->Update();

    int numBasisPts = inputRdr->GetOutput()->GetNumberOfPoints();
    int numQueryPts = numBasisPts / 2;
    std::cout << "Number of Basis Flows: " << numBasisPts/2 << "\n";
    
    // determined by macro
    std::cout << "Number of Query Pts: " << SUBSET << "\n";

    // ---------------------------------------------------------------//
    // Start timers and timing arrays
    // ---------------------------------------------------------------//
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> tri_time;
    std::chrono::duration<double> locate_time;
    std::chrono::duration<double> interp_time;
    std::vector<double> tri_times;
    std::vector<double> locate_times;
    std::vector<double> interp_times;
    std::vector<double> errors;
    // ---------------------------------------------------------------//
    // End timers and timing arrays
    // ---------------------------------------------------------------//
        
    // these hold Vertex objects
    std::vector<std::vector<Vertex>> subdomain_start_basis(64);
    std::vector<std::vector<Vertex>> subdomain_end_basis(64);
    std::vector<std::vector<Vertex>> subdomain_query_points(64);
    std::vector<Vertex> all_start_basis;

    printf("Before sorting\n");
        
    // sort into subdomains and establish bboxes
    (void)sortInputs(inputRdr, subdomain_start_basis, subdomain_end_basis, subdomain_query_points, all_start_basis, numBasisPts, doGlobal, doLocal, doNeighbor);

    // ---------------------------------------------------------------//
    // Start Interp and locate Global
    // ---------------------------------------------------------------//
    
    // destination for points and their interped end points, and a timing file
    std::ofstream outfile("results.txt");

    // keep track of locatable query points
    int unfindable = 0;

    Vertex z = Vertex(0.0,0.0,0.0);
    printf("size of Verex: %ld\n", sizeof(z));

    // For Global
    if(doGlobal){

        // Start Triang construction timer
        start = std::chrono::high_resolution_clock::now();

        // Build triangulation from ALL points
        Triangulation currentTriangulation;
        currentTriangulation.insert(all_start_basis.begin(), all_start_basis.end());

        end = std::chrono::high_resolution_clock::now();
        tri_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        tri_times.push_back(tri_time.count());
        // End Triang construction timer
        
        // check validity of triangulation
        if(!currentTriangulation.is_valid()){
            std::cerr << "TriangulationError: Global\nExiting " << argv[0] << "\n"; 
            return EXIT_FAILURE;
        } 
        
        // ---------------------------------------------------------------//
        // Start interpolation Global
        // ---------------------------------------------------------------//

        // interpolated point start
        double interpStart[3] = {0.0,0.0,0.0};

        // loop over subdomains
        for(int domId = 0; domId < 64; ++domId){
            
            // skip subdomains with no query points
            if(subdomain_query_points[domId].size() == 0) continue;

            // cell hadle grabs reference to tetrahedron vertices containing the query point
            Cell_handle cell;
            
            // Iterate over query points and then interp them 
            for(int qId = 0; qId < subdomain_query_points[domId].size(); ++qId){
                
                // Start location timer
                start = std::chrono::high_resolution_clock::now();

                // locate nearest verts via containing cell
                cell = currentTriangulation.locate(subdomain_query_points[domId][qId], cell);

                end = std::chrono::high_resolution_clock::now();
                locate_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                locate_times.push_back(locate_time.count());
                // End location timer
                
                // grab query point reference
                Vertex pt = subdomain_query_points[domId][qId];
                
                // if not valid, skip
                if(currentTriangulation.is_infinite(cell)){
                    unfindable++;
                    continue;
                }

                // grab vertices of tetrahedra containing the query point
                Vertex s0 = cell->vertex(0)->point();
                Vertex s1 = cell->vertex(1)->point();
                Vertex s2 = cell->vertex(2)->point();
                Vertex s3 = cell->vertex(3)->point();
        
                // ---------------------------------------------------------------//
                // Note: this is unecessary for this testing where start = end
                // ---------------------------------------------------------------//
                // find index of start vert in order to find corresponding end vert
                int index0 = vectorIndex(subdomain_start_basis[domId], s0);
                int index1 = vectorIndex(subdomain_start_basis[domId], s1);
                int index2 = vectorIndex(subdomain_start_basis[domId], s2);
                int index3 = vectorIndex(subdomain_start_basis[domId], s3);
                
                Vertex ev0 = subdomain_start_basis[domId][index0]; 
                Vertex ev1 = subdomain_start_basis[domId][index1]; 
                Vertex ev2 = subdomain_start_basis[domId][index2]; 
                Vertex ev3 = subdomain_start_basis[domId][index3]; 

                // declare end point for query point
                double *interpEnd = new double[3];

                // translate Vertex objects to double *
                interpStart[0] = pt.x(); 
                interpStart[1] = pt.y(); 
                interpStart[2] = pt.z(); 

                double start0[3] = {s0.x(), s0.y(), s0.z()}; 
                double start1[3] = {s1.x(), s1.y(), s1.z()}; 
                double start2[3] = {s2.x(), s2.y(), s2.z()}; 
                double start3[3] = {s3.x(), s3.y(), s3.z()}; 

                double end0[3] = {s0.x(), s0.y(), s0.z()}; 
                double end1[3] = {s1.x(), s1.y(), s1.z()};
                double end2[3] = {s2.x(), s2.y(), s2.z()};
                double end3[3] = {s3.x(), s3.y(), s3.z()};

                // DEBUG: not valid for data where start != end
                //double end0[3] = {ev0.x(), ev0.y(), ev0.z()}; 
                //double end1[3] = {ev1.x(), ev1.y(), ev1.z()}; 
                //double end2[3] = {ev2.x(), ev2.y(), ev2.z()}; 
                //double end3[3] = {ev3.x(), ev3.y(), ev3.z()}; 

                // Start interpolation timer
                start = std::chrono::high_resolution_clock::now();
                
                // interpolate
                barycentricInterp(interpEnd, interpStart, start0, start1, start2, start3, end0, end1, end2, end3);
                
                end = std::chrono::high_resolution_clock::now();
                interp_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                interp_times.push_back(interp_time.count());
                // End interpolation timing

                // calculate error (for this data set we treat end = start so error is how far we've moved)
                double error = euclidDistance3D(interpEnd, interpStart);
                errors.push_back(error); 
    
                // write start and end points to outfile
                outfile << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
                outfile << interpEnd[0] << " " << interpEnd[1] << " " << interpEnd[2] << "\n";

                // free allocated memory
                delete interpEnd;
            }
        }
        // ---------------------------------------------------------------//
        // End interpolation Global
        // ---------------------------------------------------------------//
    }

    // ---------------------------------------------------------------//
    // End Interp and locate Global
    // ---------------------------------------------------------------//

    else{
        
        // ---------------------------------------------------------------//
        // Start Interp and Locate Local and Neighbor
        // ---------------------------------------------------------------//
        double interpStart[3] = {0.0,0.0,0.0};

        // iterate over subdomains
        for(int domId = 0; domId < 64; ++domId){
            
            // need at least 4 verts for triangulation for local triang
            // we will likely have at least 4 in the neighborhood
            if(doLocal && (subdomain_start_basis[domId].size() < 4)) continue;

            // skip subdomains with no query points
            if(subdomain_query_points[domId].size() == 0) continue;
            
            // declare triangulation structure
            Triangulation currentTriangulation;

            // create cgal delaunay over current subdomain start verts only
            // Start triangulation timer
            start = std::chrono::high_resolution_clock::now();

            // load points into triangulation
            currentTriangulation.insert(subdomain_start_basis[domId].begin(), subdomain_start_basis[domId].end());

            end = std::chrono::high_resolution_clock::now();
            tri_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            tri_times.push_back(tri_time.count());
            // End triangulation timer
            
            // create cgal delaunay over current subdomain and neighbor domains
           // else if(doNeighbor){

           //     // hold neighbors
           //     std::vector<Vertex> neighborhood;
           //     buildNeighborhood(neighborhood, subdomain_start_basis, bboxes, domId);
           //     //if(true) return 0;
           //     //
           //     //// add the central subdomain
           //     //neighborhood.insert(neighborhood.begin(), subdomain_start_basis[domId].begin(), subdomain_start_basis[domId].end());

           //     //// grab all IDs for neighbor domains
           //     //std::vector<int> neighborhoodIDs = neighborDoms(subdomain_start_basis[domId][0]);

           //     //// add the neighboring subdomains
           //     ////int count = 0;
           //     //for(int i = 0; i < neighborhoodIDs.size(); ++i){

           //     //    int id = neighborhoodIDs[i];
           //     //    
           //     //    if(id != -1 && subdomain_start_basis[id].size() != 0){
           //     //            neighborhood.insert(neighborhood.begin(), subdomain_start_basis[id].begin(), subdomain_start_basis[id].end());
           //     //    }
           //     //}

           //     // Start triang timer
           //     start = std::chrono::high_resolution_clock::now();
           //     
           //     // delaunay over neighborhood
           //     currentTriangulation.insert(neighborhood.begin(), neighborhood.end());

           //     end = std::chrono::high_resolution_clock::now();
           //     tri_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
           //     tri_times.push_back(tri_time.count());
           //     // End triang timer
           // }

            // check validity of triangulation
            if(!currentTriangulation.is_valid()){
                std::cout << "Triangulation Error: domain " << domId << ": Invalid Triangulation\n"; 
                continue;
            }
           
            // ------------------------------------------------------------//
            // Start Local and Neighborhood locate and interp 
            // ------------------------------------------------------------//
            Cell_handle cell;
            for(int qId = 0; qId < subdomain_query_points[domId].size(); ++qId){

                // Start locate timer
                start = std::chrono::high_resolution_clock::now();
                
                // locate nearest verts via containing cell
                cell = currentTriangulation.locate(subdomain_query_points[domId][qId], cell);

                end = std::chrono::high_resolution_clock::now();
                locate_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                locate_times.push_back(locate_time.count());
                // End locate timer

                // grab query point
                Vertex pt = subdomain_query_points[domId][qId];
                
                // if not valid, skip
                if(currentTriangulation.is_infinite(cell)){
                    unfindable++;
                    continue;
                }


                Vertex s0 = cell->vertex(0)->point();
                Vertex s1 = cell->vertex(1)->point();
                Vertex s2 = cell->vertex(2)->point();
                Vertex s3 = cell->vertex(3)->point();

                // ---------------------------------------------------------------//
                // Note: this is unecessary for this testing where start = end
                // ---------------------------------------------------------------//
                // find index of start vert in order to find corresponding end vert
                int index0 = vectorIndex(subdomain_start_basis[domId], s0);
                int index1 = vectorIndex(subdomain_start_basis[domId], s1);
                int index2 = vectorIndex(subdomain_start_basis[domId], s2);
                int index3 = vectorIndex(subdomain_start_basis[domId], s3);
                
                // now get end vertices
                Vertex ev0 = subdomain_start_basis[domId][index0]; 
                Vertex ev1 = subdomain_start_basis[domId][index1]; 
                Vertex ev2 = subdomain_start_basis[domId][index2]; 
                Vertex ev3 = subdomain_start_basis[domId][index3]; 

                // declare end point for query point
                double *interpEnd = new double[3];

                interpStart[0] = pt.x(); 
                interpStart[1] = pt.y(); 
                interpStart[2] = pt.z(); 

                double start0[3] = {s0.x(), s0.y(), s0.z()}; 
                double start1[3] = {s1.x(), s1.y(), s1.z()}; 
                double start2[3] = {s2.x(), s2.y(), s2.z()}; 
                double start3[3] = {s3.x(), s3.y(), s3.z()}; 

                double end0[3] = {s0.x(), s0.y(), s0.z()};
                double end1[3] = {s1.x(), s1.y(), s1.z()};
                double end2[3] = {s2.x(), s2.y(), s2.z()};
                double end3[3] = {s3.x(), s3.y(), s3.z()};


                // DEBUG: not valid for data where start != end
                //double end0[3] = {ev0.x(), ev0.y(), ev0.z()}; 
                //double end1[3] = {ev1.x(), ev1.y(), ev1.z()}; 
                //double end2[3] = {ev2.x(), ev2.y(), ev2.z()}; 
                //double end3[3] = {ev3.x(), ev3.y(), ev3.z()}; 

                // Start interpolation timer
                start = std::chrono::high_resolution_clock::now();
                
                // interpolate
                barycentricInterp(interpEnd, interpStart, start0, start1, start2, start3, end0, end1, end2, end3);
                
                end = std::chrono::high_resolution_clock::now();
                interp_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                interp_times.push_back(interp_time.count());
                // End interpolation timer
                
                // calculate the error (for this data set we treat end = start so error is how far we've moved)
                double error = euclidDistance3D(interpEnd, interpStart);
                errors.push_back(error); 
                
                // write start and end points to outfile
                outfile << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
                outfile << interpEnd[0] << " " << interpEnd[1] << " " << interpEnd[2] << "\n";
                
                // free allocated memory 
                delete interpEnd;
            }
        }
    }
    // ---------------------------------------------------------------//
    // End Interp and locate Local and Neighborhood 
    // ---------------------------------------------------------------//

    // ---------------------------------------------------------------//
    // Start Print Data 
    // ---------------------------------------------------------------//

    // number of un-locatable vertices
    std::cout << "un-locatable: " << unfindable << "\n";

    // init max, min, and avg vals
    double max = tri_times[0], min = tri_times[0], avg = tri_times[0];

    std::cout << "FORMAT:\nTriangulation Time\nLocate Time(s)\nInterp Time(s)\nError(s)\n";

    std::cout << "=================================\n";
    
    // ----------------------------------------------------------- //
    // ----------------------------------------------------------- //

    if(doLocal || doNeighbor){
        avg = 0;
        for(int i = 0; i < tri_times.size(); ++i){
            if(tri_times[i] > max) max = tri_times[i];
            if(tri_times[i] < min) min = tri_times[i];
            avg += tri_times[i];
        }
        avg /= (double)tri_times.size();
    }

    std::cout << "Min: " << min << "\nMax: " << max << "\nAvg: " << avg << "\n";
    std::cout << "=================================\n";

    // ----------------------------------------------------------- //
    // ----------------------------------------------------------- //

    min = locate_times[0];
    max = locate_times[0];
    avg = 0;
    for(int i = 0; i < locate_times.size(); ++i){
        if(locate_times[i] > max) max = locate_times[i];
        if(locate_times[i] < min) min = locate_times[i];
        avg += locate_times[i];
    }
    avg /= (double)SUBSET;
    std::cout << "Min: " << min << "\nMax: " << max << "\nAvg: " << avg << "\n";
    std::cout << "=================================\n";
    
    // ----------------------------------------------------------- //
    // ----------------------------------------------------------- //

    min = interp_times[0];
    max = interp_times[0];
    avg = 0;
    for(int i = 0; i < interp_times.size(); ++i){
        if(interp_times[i] > max) max = interp_times[i];
        if(interp_times[i] < min) min = interp_times[i];
        avg += interp_times[i];
    }
    avg /= (double)SUBSET;
    std::cout << "Min: " << min << "\nMax: " << max << "\nAvg: " << avg << "\n";

    // ----------------------------------------------------------- //
    // ----------------------------------------------------------- //

    std::cout << "=================================\n";
    min = errors[0];
    max = errors[0];
    avg = 0;
    for(int i = 0; i < errors.size(); ++i){
        if(errors[i] > max) max = errors[i];
        if(errors[i] < min) min = errors[i];
        avg += errors[i];
    }
    avg /= (double)SUBSET;
    std::cout << "Min: " << min << "\nMax: " << max << "\nAvg: " << avg << "\n";
    std::cout << "=================================\n";

    // ---------------------------------------------------------------//
    // End Print Data 
    // ---------------------------------------------------------------//
    
    // ---------------------------------------------------------------//
    // Start Write Data 
    // ---------------------------------------------------------------//
//    std::string mode;
//    if(doLocal) mode = "Local";
//    if(doNeighbor) mode = "Neighbor";
//    if(doGlobal) mode = "Global";

    //std::string triString = "TriangulationConstruction_";
    //triString += mode;
    //std::ofstream triangFile(triString.c_str());

    //std::string locateString = "Location_";
    //locateString += mode;
    //std::ofstream locateFile(locateString.c_str());

    //std::string interpString = "Interp_"; 
    //interpString += mode;
    //std::ofstream interpFile(interpString.c_str());

    //std::string errorString = "Error_"; 
    //errorString += mode;
    //std::ofstream errorFile(errorString.c_str());
   
    //for(int i = 0; i < tri_times.size(); ++i){
    //    triangFile << tri_times[i] << "\n";
    //}

    //for(int i = 0; i < locate_times.size(); ++i){
    //    locateFile << locate_times[i] << "\n";
    //}

    //for(int i = 0; i < interp_times.size(); ++i){
    //    interpFile << interp_times[i] << "\n";
    //}

    //for(int i = 0; i < errors.size(); ++i){
    //    errorFile << errors[i] << "\n";
    //}

    // clean up vtk memory
    inputRdr->Delete();

    return 0;
}


