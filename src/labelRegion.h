#ifndef LABELREGION_H
#define LABELREGION_H

#include <list>

using std::list;

void labelRegion( const list<int> &region, int *ptr_seg );
void recoverLabel( const list<int> &region, int *ptr_seg );


#endif