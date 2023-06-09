# message(${CMAKE_CURRENT_SOURCE_DIR})

set(SRC
  Numeric.cpp
    fullMatrix.cpp
  BasisFactory.cpp
    FuncSpaceData.cpp
  discreteFrechetDistance.cpp
  miniBasis.cpp
    nodalBasis.cpp
	polynomialBasis.cpp
	incompleteBasis.cpp
	pyramidalBasis.cpp
	BergotBasis.cpp
	orthogonalBasis.cpp
    bezierBasis.cpp
    JacobianBasis.cpp
    CondNumBasis.cpp
    pointsGenerators.cpp
    InnerVertexPlacement.cpp
  ElementType.cpp
  GaussIntegration.cpp
    GaussQuadratureLin.cpp
    GaussQuadratureTri.cpp
    GaussQuadratureQuad.cpp
    GaussQuadratureTet.cpp
    GaussQuadratureHex.cpp
    GaussQuadraturePri.cpp
    GaussQuadraturePyr.cpp
    GaussLegendreSimplex.cpp
  GaussJacobi1D.cpp
  HilbertCurve.cpp
  robustPredicates.cpp
  decasteljau.cpp
  mathEvaluator.cpp
  Iso.cpp
  approximationError.cpp
  ConjugateGradients.cpp
)
file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)


append_gmsh_src(Numeric "${SRC};${HDR}")
