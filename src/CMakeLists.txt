#TODO INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)
INCLUDE(IJMacros.txt)


#Change PROJECT_NAME to the name of your project
PROJECT(ITKVariationalRegistration)

#The following lines are required to use Dart
#ENABLE_TESTING()
#INCLUDE(Dart)


#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
#
# The packages can be specified with a version number, for example:
#
# ITK 2.8.1
# ITK 3.2.0
#
# If no version is specified, the most recent release of the package
# will be used.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

ADD_EXECUTABLE(itkVariationalRegistration VariationalRegistrationMain.cxx) 
TARGET_LINK_LIBRARIES(itkVariationalRegistration ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkVariationalRegistration2D VariationalRegistrationMain.cxx) 
TARGET_LINK_LIBRARIES(itkVariationalRegistration2D ${ITK_LIBRARIES})
set_target_properties(itkVariationalRegistration2D PROPERTIES COMPILE_FLAGS -DUSE_2D_IMPL )

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
#ADD_TEST(CreateImage ImageCopy ${CMAKE_SOURCE_DIR}/img1.png)
#ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/img1.png ${CMAKE_SOURCE_DIR}/img1.png.base.png)

#Once you're done modifying this template, you should rename it to "CMakeLists.txt"
