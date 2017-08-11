#ifndef THRESHOLDS_H
#define THRESHOLDS_H

struct Thresholds
{
    // In  nanoseconds
    const double AP_minT;
    const double del_xtalk_minT;
    const double dir_xtalk_maxT;
    
    // In  % of pe
    const double AP;
    const double dir_xtalk;
    const double del_xtalk;

    const double time_dist;

//} default_thrs {4., 2., 2., 0.4, 1.17, 0.8, 0.4};	// H2017 Luca
} default_thrs {25., 2., 2., 0.4, 1.17, 0.8, 0.4};	// H2017 Olivier
//} default_thrs {25., 2., 2., 0.5, 1.17, 0.8, 0.4};	// H2017 Olivier
//} default_thrs {25., 2., 2., 0.6, 1.17, 0.6, 0.4};	// H2017 Olivier - PDE single threshold
//} default_thrs {4., 2., 2., 0.6, 1.2, 0.6, 0.4};
//} default_thrs {10., 5., 5., 0.5, 1.2, 0.8, 0.4};	//FBK
//} default_thrs {10., 5., 5., 0.6, 1.2, 0.6, 0.4};

#endif
