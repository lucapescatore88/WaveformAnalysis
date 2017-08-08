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

//} default_thrs {4., 2., 2., 0.4, 1.17, 0.8, 0.4};
//} default_thrs {4., 2., 2., 0.6, 1.2, 0.6, 0.4};
} default_thrs {10., 5., 5., 0.5, 1.2, 0.8, 0.4};
//} default_thrs {10., 5., 5., 0.6, 1.2, 0.6, 0.4};

#endif
