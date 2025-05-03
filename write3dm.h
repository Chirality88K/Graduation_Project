#ifndef WRITE3DM_H
#define WRITE3DM_H

#include "thirdparty/opennurbs/opennurbs.h"
#include <iostream>
#include <string>

void ChiralityWrite3dmModel(const ONX_Model* model, const wchar_t* filename);

ON_3dmObjectAttributes* Internal_CreateManagedAttributes(int layer_index, const wchar_t* name);

#endif