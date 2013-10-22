#
# gccopts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Note: For the version of system compiler supported with this release,
#       refer to the Supported and Compatible Compiler List at:
#       http://www.mathworks.com/support/compilers/current_release/
#
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building gcc MEX-files
#
# Copyright 1984-2008 The MathWorks, Inc.
# $Revision: 1.43.4.15 $  $Date: 2008/11/04 19:40:06 $
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    if [ "$ENTRYPOINT" = "mexLibrary" ]; then
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lmwservices -lut"
    else  
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat"
    fi
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS -fPIC -pthread -m32"
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS -D_FILE_OFFSET_BITS=64" 
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            CLIBS="$CLIBS -lstdc++"
#           
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -D_FILE_OFFSET_BITS=64" 
            CXXFLAGS="$CXXFLAGS -fPIC -pthread"
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
#
            # FortrankeyName: g95
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='g95'
            FFLAGS='-fexceptions'
            FFLAGS="$FFLAGS -fPIC"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexglx'
            LDFLAGS="-pthread -shared -m32 -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnxa64)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            CLIBS="$CLIBS -lstdc++"
#
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: g95
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
#
            FC='g95'
            FFLAGS='-fexceptions'
            FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexa64'
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol64)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc'
            GCC_LIBDIR=`$CC -print-file-name=libgcc_s.so | sed -e 's|libgcc_s.so||'`
            CFLAGS='-fPIC -fexceptions -m64'
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'  
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion:
            CXXDEBUGFLAGS='-g'
#
            CXX='g++'
            CXXFLAGS='-fPIC -m64'
            CXXLIBS="$MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
#
            LD="$COMPILER"
            LDEXTENSION='.mexs64'
            LDFLAGS="-shared -Wl,-M,$TMW_ROOT/extern/lib/$Arch/$MAPFILE,-R,$GCC_LIBDIR -m64"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'  
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh mac 12
#----------------------------------------------------------------------------
            ;;
        maci)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc-4.0'
            SDKROOT='/Developer/SDKs/MacOSX10.5.sdk'
            MACOSX_DEPLOYMENT_TARGET='10.5'
            ARCHS='i386'
            CFLAGS="-fno-common -no-cpp-precomp -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CFLAGS="$CFLAGS  -fexceptions"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lstdc++"
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX=g++-4.0
            CXXFLAGS="-fno-common -no-cpp-precomp -fexceptions -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: GNU Fortran
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='gfortran'
            FFLAGS='-fexceptions -fbackslash'
            FC_LIBDIR=`$FC -print-file-name=libgfortran.dylib 2>&1 | sed -n '1s/\/*libgfortran\.dylib//p'`
            FC_LIBDIR2=`$FC -print-file-name=libgfortranbegin.a 2>&1 | sed -n '1s/\/*libgfortranbegin\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lgfortran -L$FC_LIBDIR2 -lgfortranbegin"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-gdwarf-2'
#
            LD="$CC"
            LDEXTENSION='.mexmaci'
            LDFLAGS="-Wl,-twolevel_namespace -undefined error -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            LDFLAGS="$LDFLAGS -bundle -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        maci64)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc'
            SDKROOT='/Developer/SDKs/MacOSX10.6.sdk'
            MACOSX_DEPLOYMENT_TARGET='10.5'
            ARCHS='x86_64'
            CFLAGS="-fno-common -no-cpp-precomp -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CFLAGS="$CFLAGS  -fexceptions"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lstdc++"
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX=g++
			# -isysroot $SDKROOT  -> fait planter
			# CXXFLAGS='-DNDEBUG -DUNIX_MODE -DMEXMODE -fopenmp' à ajouter 
            CXXFLAGS="-fno-common -no-cpp-precomp -fexceptions -arch $ARCHS -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
			# -shared à ajouter ?
            CXXLIBS="$MLIBS -lstdc++ -shared" 
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='' # -g
  
            # FortrankeyName: GNU Fortran
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='gfortran'
            FFLAGS='-fexceptions -m64 -fbackslash'
            FC_LIBDIR=`$FC -print-file-name=libgfortran.dylib 2>&1 | sed -n '1s/\/*libgfortran\.dylib//p'`
            FC_LIBDIR2=`$FC -print-file-name=libgfortranbegin.a 2>&1 | sed -n '1s/\/*libgfortranbegin\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lgfortran -L$FC_LIBDIR2 -lgfortranbegin"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmaci64'
            LDFLAGS="-Wl,-twolevel_namespace -undefined error -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            LDFLAGS="$LDFLAGS -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -fopenmp"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
