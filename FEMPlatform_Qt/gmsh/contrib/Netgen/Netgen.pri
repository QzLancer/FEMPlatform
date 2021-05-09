INCLUDEPATH += $$PWD
INCLUDEPATH += $$PWD/nglib

win32-msvc* {
    DEFINES += _AFXDLL
}

include($$PWD/libsrc/general/general.pri)
include($$PWD/libsrc/include/include.pri)
include($$PWD/libsrc/linalg/linalg.pri)
include($$PWD/libsrc/gprim/gprim.pri)
include($$PWD/libsrc/meshing/meshing.pri)


HEADERS += \
    $$PWD/nglib/nglib.h \
    $$PWD/nglib_gmsh.h

SOURCES += \
#    $$PWD/nglib/nglib.cpp \
    $$PWD/nglib_gmsh.cpp

