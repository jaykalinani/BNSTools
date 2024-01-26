#ifndef SEPARATION_HXX
#define SEPARATION_HXX

namespace BNSTrackerX {

CCTK_DEVICE CCTK_HOST CCTK_REAL proper_distance(const cactus_glop& glop, //TODO
                        const CCTK_REAL x0, const CCTK_REAL y0, 
                        const CCTK_REAL x1, const CCTK_REAL y1, int num_points);

                        
}

#endif
 
