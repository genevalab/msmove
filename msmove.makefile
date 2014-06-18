# Compiler flags...
CPP_COMPILER = g++
C_COMPILER = gcc

# Include paths...
Debug_Include_Path=
Release_Include_Path=

# Library paths...
Debug_Library_Path=
Release_Library_Path=

# Additional libraries...
Debug_Libraries=
Release_Libraries=

# Preprocessor definitions...
Debug_Preprocessor_Definitions=-D GCC_BUILD -D _DEBUG -D _WINDOWS -D _CRT_SECURE_NO_WARNINGS 
Release_Preprocessor_Definitions=-D GCC_BUILD -D NDEBUG -D _WINDOWS -D _CRT_SECURE_NO_WARNINGS 

# Implictly linked object files...
Debug_Implicitly_Linked_Objects=
Release_Implicitly_Linked_Objects=

# Compiler flags...
Debug_Compiler_Flags=-O0 -g 
Release_Compiler_Flags=

# Builds all configurations for this project...
.PHONY: build_all_configurations
build_all_configurations: Debug Release 

# Builds the Debug configuration...
.PHONY: Debug
Debug: create_folders gccDebug/ms.o gccDebug/rand2t.o gccDebug/streec.o 
	g++ gccDebug/ms.o gccDebug/rand2t.o gccDebug/streec.o  $(Debug_Library_Path) $(Debug_Libraries) -Wl,-rpath,./ -o gccDebug/msmove

# Compiles file ms.c for the Debug configuration...
-include gccDebug/ms.d
gccDebug/ms.o: ms.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c ms.c $(Debug_Include_Path) -o gccDebug/ms.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM ms.c $(Debug_Include_Path) > gccDebug/ms.d

# Compiles file rand2t.c for the Debug configuration...
-include gccDebug/rand2t.d
gccDebug/rand2t.o: rand2t.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c rand2t.c $(Debug_Include_Path) -o gccDebug/rand2t.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM rand2t.c $(Debug_Include_Path) > gccDebug/rand2t.d

# Compiles file streec.c for the Debug configuration...
-include gccDebug/streec.d
gccDebug/streec.o: streec.c
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c streec.c $(Debug_Include_Path) -o gccDebug/streec.o
	$(C_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM streec.c $(Debug_Include_Path) > gccDebug/streec.d

# Builds the Release configuration...
.PHONY: Release
Release: create_folders gccRelease/ms.o gccRelease/rand2t.o gccRelease/streec.o 
	g++ gccRelease/ms.o gccRelease/rand2t.o gccRelease/streec.o  $(Release_Library_Path) $(Release_Libraries) -Wl,-rpath,./ -o gccRelease/msmove

# Compiles file ms.c for the Release configuration...
-include gccRelease/ms.d
gccRelease/ms.o: ms.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c ms.c $(Release_Include_Path) -o gccRelease/ms.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM ms.c $(Release_Include_Path) > gccRelease/ms.d

# Compiles file rand2t.c for the Release configuration...
-include gccRelease/rand2t.d
gccRelease/rand2t.o: rand2t.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c rand2t.c $(Release_Include_Path) -o gccRelease/rand2t.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM rand2t.c $(Release_Include_Path) > gccRelease/rand2t.d

# Compiles file streec.c for the Release configuration...
-include gccRelease/streec.d
gccRelease/streec.o: streec.c
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c streec.c $(Release_Include_Path) -o gccRelease/streec.o
	$(C_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM streec.c $(Release_Include_Path) > gccRelease/streec.d

# Creates the intermediate and output folders for each configuration...
.PHONY: create_folders
create_folders:
	mkdir -p gccDebug
	mkdir -p gccRelease

# Cleans intermediate and output files (objects, libraries, executables)...
.PHONY: clean
clean:
	rm -f gccDebug/*.o
	rm -f gccDebug/*.d
	rm -f gccDebug/*.a
	rm -f gccDebug/*.so
	rm -f gccDebug/*.dll
	rm -f gccDebug/*.exe
	rm -f gccRelease/*.o
	rm -f gccRelease/*.d
	rm -f gccRelease/*.a
	rm -f gccRelease/*.so
	rm -f gccRelease/*.dll
	rm -f gccRelease/*.exe

