#ifndef LABELREGION_H
#define LABELREGION_H

#include <list>

using std::list;

void labelRegion( const int *region, const int &len_region, int *ptr_seg );
void recoverLabel( const int *region, const int &len_region, int *ptr_seg );


#endif