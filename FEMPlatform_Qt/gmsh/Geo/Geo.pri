INCLUDEPATH += $$PWD

HEADERS += \
    $$PWD/ACISEdge.h \
    $$PWD/ACISFace.h \
    $$PWD/ACISVertex.h \
    $$PWD/affineTransformation.h \
    $$PWD/boundaryLayersData.h \
    $$PWD/Cell.h \
    $$PWD/CellComplex.h \
    $$PWD/CGNSOptions.h \
    $$PWD/CGNSUtils.h \
    $$PWD/Chain.h \
    $$PWD/ChainComplex.h \
    $$PWD/closestPoint.h \
    $$PWD/closestVertex.h \
    $$PWD/discreteEdge.h \
    $$PWD/discreteFace.h \
    $$PWD/discreteRegion.h \
    $$PWD/discreteVertex.h \
    $$PWD/ExtrudeParams.h \
    $$PWD/findLinks.h \
    $$PWD/GEdge.h \
    $$PWD/GEdgeLoop.h \
    $$PWD/GEntity.h \
    $$PWD/Geo.h \
    $$PWD/GeoDefines.h \
    $$PWD/GeoInterpolation.h \
    $$PWD/GeomMeshMatcher.h \
    $$PWD/GeoStringInterface.h \
    $$PWD/GFace.h \
    $$PWD/ghostEdge.h \
    $$PWD/ghostFace.h \
    $$PWD/ghostRegion.h \
    $$PWD/GModel.h \
    $$PWD/GModelCreateTopologyFromMesh.h \
    $$PWD/GModelIO_ACIS.h \
    $$PWD/GModelIO_GEO.h \
    $$PWD/gmshEdge.h \
    $$PWD/gmshFace.h \
    $$PWD/gmshLevelset.h \
    $$PWD/gmshRegion.h \
    $$PWD/gmshSurface.h \
    $$PWD/gmshVertex.h \
    $$PWD/GPoint.h \
    $$PWD/GRegion.h \
    $$PWD/GVertex.h \
    $$PWD/Homology.h \
    $$PWD/intersectCurveSurface.h \
    $$PWD/MEdge.h \
    $$PWD/MEdgeHash.h \
    $$PWD/MElement.h \
    $$PWD/MElementCut.h \
    $$PWD/MElementOctree.h \
    $$PWD/MFace.h \
    $$PWD/MFaceHash.h \
    $$PWD/MHexahedron.h \
    $$PWD/MLine.h \
    $$PWD/MPoint.h \
    $$PWD/MPrism.h \
    $$PWD/MPyramid.h \
    $$PWD/MQuadrangle.h \
    $$PWD/MSubElement.h \
    $$PWD/MTetrahedron.h \
    $$PWD/MTriangle.h \
    $$PWD/MTrihedron.h \
    $$PWD/MVertex.h \
    $$PWD/MVertexBoundaryLayerData.h \
    $$PWD/MVertexRTree.h \
    $$PWD/OCCAttributes.h \
    $$PWD/OCCEdge.h \
    $$PWD/OCCFace.h \
    $$PWD/OCCRegion.h \
    $$PWD/OCCVertex.h \
    $$PWD/Pair.h \
    $$PWD/partitionEdge.h \
    $$PWD/partitionFace.h \
    $$PWD/partitionRegion.h \
    $$PWD/partitionVertex.h \
    $$PWD/Range.h \
    $$PWD/SBoundingBox3d.h \
    $$PWD/SOrientedBoundingBox.h \
    $$PWD/SPoint2.h \
    $$PWD/SPoint3.h \
    $$PWD/STensor3.h \
    $$PWD/SVector3.h \
    $$PWD/xyEdge.h \
    $$PWD/xyFace.h \
    $$PWD/GModelIO_OCC.h \

SOURCES += \
    $$PWD/ACISEdge.cpp \
    $$PWD/ACISFace.cpp \
    $$PWD/ACISVertex.cpp \
    $$PWD/affineTransformation.cpp \
    $$PWD/boundaryLayersData.cpp \
    $$PWD/Cell.cpp \
    $$PWD/CellComplex.cpp \
    $$PWD/CGNSUtils.cpp \
    $$PWD/Chain.cpp \
    $$PWD/ChainComplex.cpp \
    $$PWD/closestPoint.cpp \
    $$PWD/closestVertex.cpp \
    $$PWD/discreteEdge.cpp \
    $$PWD/discreteFace.cpp \
    $$PWD/discreteRegion.cpp \
    $$PWD/discreteVertex.cpp \
    $$PWD/ExtrudeParams.cpp \
    $$PWD/findLinks.cpp \
    $$PWD/GEdge.cpp \
    $$PWD/GEdgeLoop.cpp \
    $$PWD/GEntity.cpp \
    $$PWD/Geo.cpp \
    $$PWD/GeoInterpolation.cpp \
    $$PWD/GeomMeshMatcher.cpp \
    $$PWD/GeoStringInterface.cpp \
    $$PWD/GFace.cpp \
    $$PWD/GModel.cpp \
    $$PWD/GModelCreateTopologyFromMesh.cpp \
    $$PWD/GModelIO_ACIS.cpp \
    $$PWD/GModelIO_ACTRAN.cpp \
    $$PWD/GModelIO_BDF.cpp \
    $$PWD/GModelIO_CELUM.cpp \
    $$PWD/GModelIO_CGNS.cpp \
    $$PWD/GModelIO_DIFF.cpp \
    $$PWD/GModelIO_GEO.cpp \
    $$PWD/GModelIO_GEOM.cpp \
    $$PWD/GModelIO_INP.cpp \
    $$PWD/GModelIO_IR3.cpp \
    $$PWD/GModelIO_KEY.cpp \
    $$PWD/GModelIO_MAIL.cpp \
    $$PWD/GModelIO_MATLAB.cpp \
    $$PWD/GModelIO_MED.cpp \
    $$PWD/GModelIO_MESH.cpp \
    $$PWD/GModelIO_MSH.cpp \
    $$PWD/GModelIO_MSH2.cpp \
    $$PWD/GModelIO_MSH3.cpp \
    $$PWD/GModelIO_MSH4.cpp \
    $$PWD/GModelIO_NEU.cpp \
    $$PWD/GModelIO_P3D.cpp \
    $$PWD/GModelIO_PLY.cpp \
    $$PWD/GModelIO_POS.cpp \
    $$PWD/GModelIO_SAMCEF.cpp \
    $$PWD/GModelIO_STL.cpp \
    $$PWD/GModelIO_SU2.cpp \
    $$PWD/GModelIO_TOCHNOG.cpp \
    $$PWD/GModelIO_UNV.cpp \
    $$PWD/GModelIO_VRML.cpp \
    $$PWD/GModelIO_VTK.cpp \
    $$PWD/GModelVertexArrays.cpp \
    $$PWD/gmshEdge.cpp \
    $$PWD/gmshFace.cpp \
    $$PWD/gmshLevelset.cpp \
    $$PWD/gmshRegion.cpp \
    $$PWD/gmshSurface.cpp \
    $$PWD/gmshVertex.cpp \
    $$PWD/GRegion.cpp \
    $$PWD/GVertex.cpp \
    $$PWD/Homology.cpp \
    $$PWD/intersectCurveSurface.cpp \
    $$PWD/MEdge.cpp \
    $$PWD/MElement.cpp \
    $$PWD/MElementCut.cpp \
    $$PWD/MElementOctree.cpp \
    $$PWD/MFace.cpp \
    $$PWD/MHexahedron.cpp \
    $$PWD/MLine.cpp \
    $$PWD/MPrism.cpp \
    $$PWD/MPyramid.cpp \
    $$PWD/MQuadrangle.cpp \
    $$PWD/MSubElement.cpp \
    $$PWD/MTetrahedron.cpp \
    $$PWD/MTriangle.cpp \
    $$PWD/MTrihedron.cpp \
    $$PWD/MVertex.cpp \
    $$PWD/MVertexBoundaryLayerData.cpp \
    $$PWD/OCCEdge.cpp \
    $$PWD/OCCFace.cpp \
    $$PWD/OCCRegion.cpp \
    $$PWD/OCCVertex.cpp \
    $$PWD/SOrientedBoundingBox.cpp \
    $$PWD/STensor3.cpp \
    $$PWD/GModelIO_OCC.cpp \