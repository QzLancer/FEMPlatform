// Copyright 2020 Poofee (https://github.com/Poofee)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------
/*****************************************************************************
 *                                                                           *
 *  File:    meshGRegiontetgen.cpp                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Poofee                                                          *
 *  Email:   poofee@qq.com                                                   *
 *  Address: https://github.com/Poofee                                       *
 *  Date:    2020年11月04日                                                   *
 *                                                                           *
 *****************************************************************************/

#include <stdlib.h>
#include <vector>
#include "GmshConfig.h"
#include "GmshMessage.h"
#include "meshGRegion.h"
#include "GModel.h"
#include "GRegion.h"
#include "GFace.h"
#include "MTriangle.h"
#include "MTetrahedron.h"
#include "ExtrudeParams.h"
#include "Context.h"

#if defined(HAVE_TETGEN)


#endif

void meshGRegiontetgen(GRegion *gr)
{
#if !defined(HAVE_TETGEN)
  Msg::Error("Requires Tetgen");
#else
  // sanity check for frontal algo
  std::vector<GFace *> faces = gr->faces();
  for(std::vector<GFace *>::iterator it = faces.begin(); it != faces.end(); it++) {
    if((*it)->quadrangles.size()) {
      Msg::Error("Cannot use frontal 3D algorithm with quadrangles on boundary");
      return;
    }
  }

  Msg::Info("Meshing volume %d (Frontal)", gr->tag());
  
#endif
}

